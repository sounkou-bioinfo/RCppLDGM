#include <Rcpp.h>
#include <hdf5.h>
#include <hdf5lib.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace {

constexpr H5Z_filter_t H5Z_FILTER_LZF_RCPPLDGM = 32000;

class H5Handle {
 public:
  using Closer = herr_t (*)(hid_t);

  H5Handle() = default;
  H5Handle(hid_t id, Closer closer) : id_(id), closer_(closer) {}
  H5Handle(const H5Handle&) = delete;
  H5Handle& operator=(const H5Handle&) = delete;

  H5Handle(H5Handle&& other) noexcept : id_(other.id_), closer_(other.closer_) {
    other.id_ = -1;
    other.closer_ = nullptr;
  }

  H5Handle& operator=(H5Handle&& other) noexcept {
    if (this != &other) {
      reset();
      id_ = other.id_;
      closer_ = other.closer_;
      other.id_ = -1;
      other.closer_ = nullptr;
    }
    return *this;
  }

  ~H5Handle() { reset(); }

  hid_t get() const { return id_; }
  operator hid_t() const { return id_; }

  void reset() {
    if (id_ >= 0 && closer_ != nullptr) {
      closer_(id_);
    }
    id_ = -1;
    closer_ = nullptr;
  }

 private:
  hid_t id_ = -1;
  Closer closer_ = nullptr;
};

[[noreturn]] void h5_fail(const std::string& what) {
  Rcpp::stop("HDF5 error while " + what);
}

void check_status(herr_t status, const std::string& what) {
  if (status < 0) {
    h5_fail(what);
  }
}

hid_t check_id(hid_t id, const std::string& what) {
  if (id < 0) {
    h5_fail(what);
  }
  return id;
}

void ensure_hdf5_ready() {
  check_status(H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr), "disabling HDF5 error printing");
  static bool filters_registered = false;
  if (!filters_registered) {
    check_status(hdf5lib_register_all_filters(), "registering hdf5lib compression filters");
    filters_registered = true;
  }
}

bool link_exists(hid_t loc, const std::string& name) {
  htri_t exists = H5Lexists(loc, name.c_str(), H5P_DEFAULT);
  if (exists < 0) {
    h5_fail("checking link '" + name + "'");
  }
  return exists > 0;
}

H5Handle create_group(hid_t loc, const std::string& name) {
  return H5Handle(check_id(H5Gcreate2(loc, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                           "creating group '" + name + "'"),
                  H5Gclose);
}

H5Handle open_group(hid_t loc, const std::string& name) {
  return H5Handle(check_id(H5Gopen2(loc, name.c_str(), H5P_DEFAULT),
                           "opening group '" + name + "'"),
                  H5Gclose);
}

H5Handle require_group(hid_t loc, const std::string& name) {
  if (link_exists(loc, name)) {
    return open_group(loc, name);
  }
  return create_group(loc, name);
}

H5Handle make_vlen_string_type() {
  H5Handle type(check_id(H5Tcopy(H5T_C_S1), "copying HDF5 string type"), H5Tclose);
  check_status(H5Tset_size(type.get(), H5T_VARIABLE), "setting variable string length");
  check_status(H5Tset_cset(type.get(), H5T_CSET_UTF8), "setting UTF-8 string encoding");
  return type;
}

H5Handle make_scalar_space() {
  return H5Handle(check_id(H5Screate(H5S_SCALAR), "creating scalar dataspace"), H5Sclose);
}

H5Handle make_1d_space(hsize_t n) {
  return H5Handle(check_id(H5Screate_simple(1, &n, nullptr), "creating one-dimensional dataspace"), H5Sclose);
}

H5Handle make_dataset_plist(hsize_t n, const std::string& compression, int chunk_size) {
  H5Handle plist(check_id(H5Pcreate(H5P_DATASET_CREATE), "creating dataset property list"), H5Pclose);
  hsize_t chunk = static_cast<hsize_t>(std::max(1, chunk_size));
  if (n > 0) {
    chunk = std::min(chunk, n);
  }
  check_status(H5Pset_chunk(plist.get(), 1, &chunk), "setting HDF5 chunk size");

  if (compression == "lzf") {
    htri_t available = H5Zfilter_avail(H5Z_FILTER_LZF_RCPPLDGM);
    if (available <= 0) {
      h5_fail("checking availability of the hdf5lib LZF filter");
    }
    check_status(H5Pset_filter(plist.get(), H5Z_FILTER_LZF_RCPPLDGM, H5Z_FLAG_MANDATORY, 0, nullptr),
                 "enabling LZF compression");
  } else if (compression == "gzip") {
    check_status(H5Pset_deflate(plist.get(), 4), "enabling gzip compression");
  } else if (compression != "none") {
    Rcpp::stop("unsupported HDF5 compression '%s'", compression);
  }
  return plist;
}

std::string scalar_string(Rcpp::String value, const std::string& name) {
  if (value == NA_STRING) {
    Rcpp::stop("`%s` must not contain missing strings", name);
  }
  return std::string(value.get_cstring());
}

std::vector<std::string> character_vector_to_strings(Rcpp::CharacterVector values,
                                                     const std::string& name) {
  std::vector<std::string> out(values.size());
  for (R_xlen_t i = 0; i < values.size(); ++i) {
    out[static_cast<size_t>(i)] = scalar_string(values[i], name);
  }
  return out;
}

std::vector<const char*> string_ptrs(const std::vector<std::string>& values) {
  std::vector<const char*> ptrs(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    ptrs[i] = values[i].c_str();
  }
  return ptrs;
}

std::vector<long long> numeric_to_int64(SEXP values, const std::string& name) {
  R_xlen_t n = Rf_xlength(values);
  std::vector<long long> out(static_cast<size_t>(n));
  if (Rf_isInteger(values)) {
    Rcpp::IntegerVector x(values);
    for (R_xlen_t i = 0; i < n; ++i) {
      if (x[i] == NA_INTEGER) {
        Rcpp::stop("`%s` must not contain missing values", name);
      }
      out[static_cast<size_t>(i)] = static_cast<long long>(x[i]);
    }
  } else if (Rf_isReal(values)) {
    Rcpp::NumericVector x(values);
    const double min_ll = static_cast<double>(std::numeric_limits<long long>::min());
    const double max_ll = static_cast<double>(std::numeric_limits<long long>::max());
    for (R_xlen_t i = 0; i < n; ++i) {
      double value = x[i];
      if (!R_finite(value) || std::floor(value) != value || value < min_ll || value > max_ll) {
        Rcpp::stop("`%s` must contain finite integer-valued data", name);
      }
      out[static_cast<size_t>(i)] = static_cast<long long>(value);
    }
  } else {
    Rcpp::stop("`%s` must be integer or numeric", name);
  }
  return out;
}

void write_attr_string(hid_t loc, const std::string& name, const std::string& value) {
  if (H5Aexists(loc, name.c_str()) > 0) {
    check_status(H5Adelete(loc, name.c_str()), "deleting existing attribute '" + name + "'");
  }
  H5Handle type = make_vlen_string_type();
  H5Handle space = make_scalar_space();
  H5Handle attr(check_id(H5Acreate2(loc, name.c_str(), type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT),
                         "creating attribute '" + name + "'"),
                H5Aclose);
  const char* ptr = value.c_str();
  check_status(H5Awrite(attr.get(), type.get(), &ptr), "writing attribute '" + name + "'");
}

void write_attr_strings(hid_t loc, const std::string& name, const std::vector<std::string>& values) {
  if (H5Aexists(loc, name.c_str()) > 0) {
    check_status(H5Adelete(loc, name.c_str()), "deleting existing attribute '" + name + "'");
  }
  hsize_t n = static_cast<hsize_t>(values.size());
  H5Handle type = make_vlen_string_type();
  H5Handle space(check_id(H5Screate_simple(1, &n, nullptr), "creating string-vector attribute dataspace"),
                 H5Sclose);
  H5Handle attr(check_id(H5Acreate2(loc, name.c_str(), type.get(), space.get(), H5P_DEFAULT, H5P_DEFAULT),
                         "creating attribute '" + name + "'"),
                H5Aclose);
  std::vector<const char*> ptrs = string_ptrs(values);
  check_status(H5Awrite(attr.get(), type.get(), ptrs.data()), "writing attribute '" + name + "'");
}

void write_dataset_int64(hid_t loc,
                         const std::string& name,
                         const std::vector<long long>& values,
                         const std::string& compression,
                         int chunk_size) {
  hsize_t n = static_cast<hsize_t>(values.size());
  H5Handle space = make_1d_space(n);
  H5Handle plist = make_dataset_plist(n, compression, chunk_size);
  H5Handle dataset(check_id(H5Dcreate2(loc, name.c_str(), H5T_NATIVE_LLONG, space.get(), H5P_DEFAULT,
                                       plist.get(), H5P_DEFAULT),
                            "creating dataset '" + name + "'"),
                   H5Dclose);
  check_status(H5Dwrite(dataset.get(), H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.data()),
               "writing dataset '" + name + "'");
}

void write_dataset_double(hid_t loc,
                          const std::string& name,
                          Rcpp::NumericVector values,
                          const std::string& compression,
                          int chunk_size) {
  hsize_t n = static_cast<hsize_t>(values.size());
  H5Handle space = make_1d_space(n);
  H5Handle plist = make_dataset_plist(n, compression, chunk_size);
  H5Handle dataset(check_id(H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, space.get(), H5P_DEFAULT,
                                       plist.get(), H5P_DEFAULT),
                            "creating dataset '" + name + "'"),
                   H5Dclose);
  check_status(H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, values.begin()),
               "writing dataset '" + name + "'");
}

void write_dataset_double_matrix(hid_t loc,
                                 const std::string& name,
                                 Rcpp::NumericMatrix values,
                                 const std::string& compression,
                                 int chunk_size) {
  hsize_t dims[2] = {static_cast<hsize_t>(values.nrow()), static_cast<hsize_t>(values.ncol())};
  H5Handle space(check_id(H5Screate_simple(2, dims, nullptr), "creating two-dimensional dataspace"), H5Sclose);
  H5Handle plist(check_id(H5Pcreate(H5P_DATASET_CREATE), "creating matrix dataset property list"), H5Pclose);
  hsize_t chunk[2] = {std::max<hsize_t>(1, std::min<hsize_t>(dims[0], static_cast<hsize_t>(std::max(1, chunk_size)))),
                      std::max<hsize_t>(1, dims[1])};
  check_status(H5Pset_chunk(plist.get(), 2, chunk), "setting HDF5 matrix chunk size");
  if (compression == "lzf") {
    htri_t available = H5Zfilter_avail(H5Z_FILTER_LZF_RCPPLDGM);
    if (available <= 0) {
      h5_fail("checking availability of the hdf5lib LZF filter");
    }
    check_status(H5Pset_filter(plist.get(), H5Z_FILTER_LZF_RCPPLDGM, H5Z_FLAG_MANDATORY, 0, nullptr),
                 "enabling LZF compression");
  } else if (compression == "gzip") {
    check_status(H5Pset_deflate(plist.get(), 4), "enabling gzip compression");
  } else if (compression != "none") {
    Rcpp::stop("unsupported HDF5 compression '%s'", compression);
  }
  H5Handle dataset(check_id(H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, space.get(), H5P_DEFAULT,
                                       plist.get(), H5P_DEFAULT),
                            "creating dataset '" + name + "'"),
                   H5Dclose);
  std::vector<double> buffer(static_cast<size_t>(values.nrow()) * static_cast<size_t>(values.ncol()));
  for (int i = 0; i < values.nrow(); ++i) {
    for (int j = 0; j < values.ncol(); ++j) {
      buffer[static_cast<size_t>(i) * static_cast<size_t>(values.ncol()) + static_cast<size_t>(j)] = values(i, j);
    }
  }
  check_status(H5Dwrite(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()),
               "writing dataset '" + name + "'");
}

void write_dataset_strings(hid_t loc,
                           const std::string& name,
                           const std::vector<std::string>& values,
                           const std::string& compression,
                           int chunk_size) {
  hsize_t n = static_cast<hsize_t>(values.size());
  H5Handle type = make_vlen_string_type();
  H5Handle space = make_1d_space(n);
  H5Handle plist = make_dataset_plist(n, compression, chunk_size);
  H5Handle dataset(check_id(H5Dcreate2(loc, name.c_str(), type.get(), space.get(), H5P_DEFAULT,
                                       plist.get(), H5P_DEFAULT),
                            "creating dataset '" + name + "'"),
                   H5Dclose);
  std::vector<const char*> ptrs = string_ptrs(values);
  check_status(H5Dwrite(dataset.get(), type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.data()),
               "writing dataset '" + name + "'");
}

H5Handle open_file_for_write(const std::string& filename, bool overwrite) {
  if (overwrite) {
    return H5Handle(check_id(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
                             "creating HDF5 file '" + filename + "'"),
                    H5Fclose);
  }

  FILE* probe = std::fopen(filename.c_str(), "rb");
  if (probe != nullptr) {
    std::fclose(probe);
    return H5Handle(check_id(H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT),
                             "opening HDF5 file '" + filename + "' for update"),
                    H5Fclose);
  }
  return H5Handle(check_id(H5Fcreate(filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT),
                           "creating HDF5 file '" + filename + "'"),
                  H5Fclose);
}

void create_variant_group_if_needed(hid_t file,
                                    Rcpp::DataFrame variant_data,
                                    Rcpp::IntegerVector jackknife_blocks,
                                    const std::string& source,
                                    const std::string& compression,
                                    int chunk_size) {
  if (link_exists(file, "row_data")) {
    require_group(file, "traits");
    require_group(file, "groups");
    return;
  }

  Rcpp::CharacterVector rsid = variant_data["RSID"];
  SEXP chr = variant_data["CHR"];
  SEXP pos = variant_data["POS"];

  write_attr_string(file, "metadata", "");
  write_attr_string(file, "data_type", "variant");
  write_attr_strings(file, "keys", std::vector<std::string>{"RSID", "POS", "CHR"});
  write_attr_string(file, "source", source);

  require_group(file, "traits");
  H5Handle row_data = create_group(file, "row_data");
  require_group(file, "groups");

  write_dataset_int64(row_data.get(), "CHR", numeric_to_int64(chr, "variant_data$CHR"), compression, chunk_size);
  write_dataset_int64(row_data.get(), "POS", numeric_to_int64(pos, "variant_data$POS"), compression, chunk_size);
  write_dataset_strings(row_data.get(), "RSID", character_vector_to_strings(rsid, "variant_data$RSID"), compression,
                        chunk_size);
  write_dataset_int64(row_data.get(), "jackknife_blocks", numeric_to_int64(jackknife_blocks, "jackknife_blocks"),
                      compression, chunk_size);
}

std::vector<hsize_t> dataset_dims(hid_t dataset) {
  H5Handle space(check_id(H5Dget_space(dataset), "opening dataset dataspace"), H5Sclose);
  int rank = H5Sget_simple_extent_ndims(space.get());
  if (rank < 0) {
    h5_fail("reading dataset rank");
  }
  std::vector<hsize_t> dims(static_cast<size_t>(rank));
  if (rank > 0) {
    check_status(H5Sget_simple_extent_dims(space.get(), dims.data(), nullptr), "reading dataset dimensions");
  }
  return dims;
}

hsize_t flat_dataset_length(hid_t dataset) {
  std::vector<hsize_t> dims = dataset_dims(dataset);
  if (dims.empty()) {
    return 1;
  }
  hsize_t n = 1;
  for (hsize_t dim : dims) {
    n *= dim;
  }
  return n;
}

Rcpp::NumericVector read_dataset_int64_as_numeric(hid_t loc, const std::string& name) {
  H5Handle dataset(check_id(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "opening dataset '" + name + "'"), H5Dclose);
  hsize_t n = flat_dataset_length(dataset.get());
  std::vector<long long> buffer(static_cast<size_t>(n));
  check_status(H5Dread(dataset.get(), H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()),
               "reading dataset '" + name + "'");
  Rcpp::NumericVector out(static_cast<R_xlen_t>(n));
  for (hsize_t i = 0; i < n; ++i) {
    out[static_cast<R_xlen_t>(i)] = static_cast<double>(buffer[static_cast<size_t>(i)]);
  }
  return out;
}

Rcpp::NumericVector read_dataset_double(hid_t loc, const std::string& name) {
  H5Handle dataset(check_id(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "opening dataset '" + name + "'"), H5Dclose);
  hsize_t n = flat_dataset_length(dataset.get());
  Rcpp::NumericVector out(static_cast<R_xlen_t>(n));
  check_status(H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.begin()),
               "reading dataset '" + name + "'");
  return out;
}

Rcpp::NumericMatrix read_dataset_double_matrix(hid_t loc, const std::string& name) {
  H5Handle dataset(check_id(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "opening dataset '" + name + "'"), H5Dclose);
  std::vector<hsize_t> dims = dataset_dims(dataset.get());
  if (dims.size() != 2) {
    Rcpp::stop("dataset '%s' must be two-dimensional", name);
  }
  Rcpp::NumericMatrix out(static_cast<int>(dims[0]), static_cast<int>(dims[1]));
  std::vector<double> buffer(static_cast<size_t>(dims[0]) * static_cast<size_t>(dims[1]));
  check_status(H5Dread(dataset.get(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()),
               "reading dataset '" + name + "'");
  for (hsize_t i = 0; i < dims[0]; ++i) {
    for (hsize_t j = 0; j < dims[1]; ++j) {
      out(static_cast<int>(i), static_cast<int>(j)) = buffer[static_cast<size_t>(i) * static_cast<size_t>(dims[1]) + static_cast<size_t>(j)];
    }
  }
  return out;
}

Rcpp::CharacterVector read_dataset_strings(hid_t loc, const std::string& name) {
  H5Handle dataset(check_id(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "opening dataset '" + name + "'"), H5Dclose);
  hsize_t n = flat_dataset_length(dataset.get());
  H5Handle type(check_id(H5Dget_type(dataset.get()), "opening datatype for dataset '" + name + "'"), H5Tclose);
  H5Handle space(check_id(H5Dget_space(dataset.get()), "opening dataspace for dataset '" + name + "'"), H5Sclose);
  Rcpp::CharacterVector out(static_cast<R_xlen_t>(n));

  if (H5Tis_variable_str(type.get()) > 0) {
    std::vector<char*> buffer(static_cast<size_t>(n), nullptr);
    check_status(H5Dread(dataset.get(), type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()),
                 "reading string dataset '" + name + "'");
    for (hsize_t i = 0; i < n; ++i) {
      char* ptr = buffer[static_cast<size_t>(i)];
      out[static_cast<R_xlen_t>(i)] = ptr == nullptr ? "" : ptr;
    }
    check_status(H5Dvlen_reclaim(type.get(), space.get(), H5P_DEFAULT, buffer.data()),
                 "reclaiming variable-length string memory");
  } else {
    size_t width = H5Tget_size(type.get());
    if (width == 0) {
      h5_fail("reading fixed-width string size for dataset '" + name + "'");
    }
    std::vector<char> buffer(static_cast<size_t>(n) * width, '\0');
    check_status(H5Dread(dataset.get(), type.get(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer.data()),
                 "reading fixed-width string dataset '" + name + "'");
    for (hsize_t i = 0; i < n; ++i) {
      const char* start = buffer.data() + static_cast<size_t>(i) * width;
      size_t len = 0;
      while (len < width && start[len] != '\0') {
        ++len;
      }
      out[static_cast<R_xlen_t>(i)] = std::string(start, len);
    }
  }
  return out;
}

Rcpp::CharacterVector list_group_names(hid_t group) {
  H5G_info_t info;
  check_status(H5Gget_info(group, &info), "reading group link count");
  Rcpp::CharacterVector out(static_cast<R_xlen_t>(info.nlinks));
  for (hsize_t i = 0; i < info.nlinks; ++i) {
    ssize_t needed = H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, i, nullptr, 0, H5P_DEFAULT);
    if (needed < 0) {
      h5_fail("reading group link name length");
    }
    std::vector<char> buffer(static_cast<size_t>(needed) + 1, '\0');
    ssize_t written = H5Lget_name_by_idx(group, ".", H5_INDEX_NAME, H5_ITER_INC, i, buffer.data(), buffer.size(),
                                         H5P_DEFAULT);
    if (written < 0) {
      h5_fail("reading group link name");
    }
    out[static_cast<R_xlen_t>(i)] = std::string(buffer.data(), static_cast<size_t>(written));
  }
  return out;
}

}  // namespace

// [[Rcpp::export(name = "RC_hdf5_filter_info")]]
Rcpp::List hdf5_filter_info_cpp() {
  ensure_hdf5_ready();
  return Rcpp::List::create(
      Rcpp::Named("lzf") = (H5Zfilter_avail(H5Z_FILTER_LZF_RCPPLDGM) > 0),
      Rcpp::Named("gzip") = (H5Zfilter_avail(H5Z_FILTER_DEFLATE) > 0));
}

// [[Rcpp::export(name = "RC_write_graphld_score_hdf5")]]
Rcpp::List write_graphld_score_hdf5_cpp(const std::string& filename,
                                        Rcpp::DataFrame variant_data,
                                        Rcpp::NumericVector gradient,
                                        SEXP hessian,
                                        const std::string& trait_name,
                                        Rcpp::IntegerVector jackknife_blocks,
                                        bool overwrite,
                                        const std::string& source,
                                        const std::string& compression,
                                        int chunk_size,
                                        SEXP parameters,
                                        SEXP jackknife_parameters) {
  ensure_hdf5_ready();
  if (gradient.size() == 0) {
    Rcpp::stop("`gradient` must contain at least one value");
  }
  if (trait_name.empty()) {
    Rcpp::stop("`trait_name` must not be empty");
  }
  if (trait_name.find('/') != std::string::npos) {
    Rcpp::stop("`trait_name` must not contain '/'");
  }
  if (jackknife_blocks.size() != gradient.size()) {
    Rcpp::stop("`jackknife_blocks` length must equal `gradient` length");
  }
  if (chunk_size < 1) {
    Rcpp::stop("`chunk_size` must be positive");
  }

  H5Handle file = open_file_for_write(filename, overwrite);
  create_variant_group_if_needed(file.get(), variant_data, jackknife_blocks, source, compression, chunk_size);
  H5Handle row_data = open_group(file.get(), "row_data");
  H5Handle rsid_dataset(check_id(H5Dopen2(row_data.get(), "RSID", H5P_DEFAULT),
                                 "opening dataset 'row_data/RSID'"),
                        H5Dclose);
  if (flat_dataset_length(rsid_dataset.get()) != static_cast<hsize_t>(gradient.size())) {
    Rcpp::stop("existing HDF5 row_data length must equal `gradient` length");
  }

  H5Handle traits = require_group(file.get(), "traits");
  if (link_exists(traits.get(), trait_name)) {
    Rcpp::stop("the HDF5 group 'traits/%s' already exists", trait_name);
  }
  H5Handle trait_group = create_group(traits.get(), trait_name);
  write_dataset_double(trait_group.get(), "gradient", gradient, compression, chunk_size);

  bool wrote_hessian = !Rf_isNull(hessian);
  if (wrote_hessian) {
    Rcpp::NumericVector hessian_values(hessian);
    if (hessian_values.size() != gradient.size()) {
      Rcpp::stop("`hessian` length must equal `gradient` length");
    }
    write_dataset_double(trait_group.get(), "hessian", hessian_values, compression, chunk_size);
  }

  bool wrote_parameters = !Rf_isNull(parameters);
  if (wrote_parameters) {
    if (Rf_isNull(jackknife_parameters)) {
      Rcpp::stop("`jackknife_parameters` must be supplied when `parameters` is supplied");
    }
    Rcpp::NumericVector parameter_values(parameters);
    Rcpp::NumericMatrix jackknife_values(jackknife_parameters);
    if (jackknife_values.ncol() != parameter_values.size()) {
      Rcpp::stop("`jackknife_parameters` column count must equal `parameters` length");
    }
    H5Handle parameter_group = create_group(trait_group.get(), "parameters");
    write_dataset_double(parameter_group.get(), "parameters", parameter_values, compression, chunk_size);
    write_dataset_double_matrix(parameter_group.get(), "jackknife_parameters", jackknife_values, compression,
                                chunk_size);
  }
  check_status(H5Fflush(file.get(), H5F_SCOPE_GLOBAL), "flushing HDF5 file");

  return Rcpp::List::create(Rcpp::Named("file") = filename,
                            Rcpp::Named("trait_name") = trait_name,
                            Rcpp::Named("n_variants") = gradient.size(),
                            Rcpp::Named("compression") = compression,
                            Rcpp::Named("hessian") = wrote_hessian,
                            Rcpp::Named("parameters") = wrote_parameters);
}

// [[Rcpp::export(name = "RC_read_graphld_score_hdf5")]]
Rcpp::List read_graphld_score_hdf5_cpp(const std::string& filename, const std::string& trait_name) {
  ensure_hdf5_ready();
  H5Handle file(check_id(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
                         "opening HDF5 file '" + filename + "' for reading"),
                H5Fclose);
  H5Handle row_data = open_group(file.get(), "row_data");
  H5Handle traits = open_group(file.get(), "traits");

  Rcpp::DataFrame variant_data = Rcpp::DataFrame::create(
      Rcpp::Named("CHR") = read_dataset_int64_as_numeric(row_data.get(), "CHR"),
      Rcpp::Named("POS") = read_dataset_int64_as_numeric(row_data.get(), "POS"),
      Rcpp::Named("RSID") = read_dataset_strings(row_data.get(), "RSID"),
      Rcpp::Named("jackknife_blocks") = read_dataset_int64_as_numeric(row_data.get(), "jackknife_blocks"),
      Rcpp::_["stringsAsFactors"] = false);

  Rcpp::CharacterVector trait_names = list_group_names(traits.get());
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("variant_data") = variant_data,
                                      Rcpp::Named("trait_names") = trait_names);
  if (!trait_name.empty()) {
    if (!link_exists(traits.get(), trait_name)) {
      Rcpp::stop("the HDF5 group 'traits/%s' does not exist", trait_name);
    }
    H5Handle trait_group = open_group(traits.get(), trait_name);
    out["trait_name"] = trait_name;
    out["gradient"] = read_dataset_double(trait_group.get(), "gradient");
    if (link_exists(trait_group.get(), "hessian")) {
      out["hessian"] = read_dataset_double(trait_group.get(), "hessian");
    }
    if (link_exists(trait_group.get(), "parameters")) {
      H5Handle parameter_group = open_group(trait_group.get(), "parameters");
      out["parameters"] = read_dataset_double(parameter_group.get(), "parameters");
      out["jackknife_parameters"] = read_dataset_double_matrix(parameter_group.get(), "jackknife_parameters");
    }
  }
  return out;
}

// [[Rcpp::export(name = "RC_write_graphld_surrogate_hdf5")]]
Rcpp::List write_graphld_surrogate_hdf5_cpp(const std::string& filename,
                                            const std::string& block_name,
                                            SEXP surrogate_map,
                                            bool overwrite,
                                            const std::string& compression,
                                            int chunk_size) {
  ensure_hdf5_ready();
  if (block_name.empty()) {
    Rcpp::stop("`block_name` must not be empty");
  }
  if (block_name.find('/') != std::string::npos) {
    Rcpp::stop("`block_name` must not contain '/'");
  }
  if (chunk_size < 1) {
    Rcpp::stop("`chunk_size` must be positive");
  }
  std::vector<long long> values = numeric_to_int64(surrogate_map, "surrogate_map");
  if (values.empty()) {
    Rcpp::stop("`surrogate_map` must contain at least one value");
  }

  H5Handle file = open_file_for_write(filename, overwrite);
  if (link_exists(file.get(), block_name)) {
    Rcpp::stop("the HDF5 dataset '%s' already exists", block_name);
  }
  write_dataset_int64(file.get(), block_name, values, compression, chunk_size);
  check_status(H5Fflush(file.get(), H5F_SCOPE_GLOBAL), "flushing HDF5 file");
  return Rcpp::List::create(Rcpp::Named("file") = filename,
                            Rcpp::Named("block_name") = block_name,
                            Rcpp::Named("n") = static_cast<int>(values.size()),
                            Rcpp::Named("compression") = compression);
}

// [[Rcpp::export(name = "RC_read_graphld_surrogate_hdf5")]]
Rcpp::NumericVector read_graphld_surrogate_hdf5_cpp(const std::string& filename,
                                                    const std::string& block_name) {
  ensure_hdf5_ready();
  if (block_name.empty()) {
    Rcpp::stop("`block_name` must not be empty");
  }
  H5Handle file(check_id(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
                         "opening HDF5 file '" + filename + "' for reading"),
                H5Fclose);
  if (!link_exists(file.get(), block_name)) {
    Rcpp::stop("the HDF5 dataset '%s' does not exist", block_name);
  }
  return read_dataset_int64_as_numeric(file.get(), block_name);
}
