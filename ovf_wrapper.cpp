#include <iostream>
#include "ovf_wrapper.hpp"
#include "utils.h"

namespace ovf {

ovf_segment *Segment::operator->() {
    return ptr.get();
}

const ovf_segment *Segment::operator->() const {
    return ptr.get();
}

ovf_file *File::operator->() const {
    return ptr;
}

File File::open_file(const char *path) {
    auto ptr = ovf_open(path);
    if (!ptr->found) {
        throw std::runtime_error(join_string("ovf file `", path, "` not found"));
    }
    return File(ptr);
}

File File::create_file(const char *path) {
    { std::ofstream _create_file(path); }
    return open_file(path);
}

File::~File() {
    handle_err("file destructor", ovf_close);
    delete ptr;
}

arma::mat File::read(int idx) {
    Segment seg;
    handle_err("read segment header", ovf_read_segment_header, idx, seg.ptr.get());
    // have to flip rows <-> cols, then transpose
    mat ret = seg.is_rectangular()
              ? mat(seg->valuedim, seg->N)
              : mat(seg->valuedim + 3, seg->pointcount);
    handle_err("read segment data", ovf_read_segment_data_8, idx, seg.ptr.get(), ret.memptr());
    if (segments.size() < idx + 1) {
        segments.resize(idx + 1);
    }
    segments[idx] = std::move(seg);
    return ret.t();
}

std::vector<arma::mat> File::read_all() {
    std::vector<arma::mat> ret;
    for (int i = 0; i < ptr->n_segments; ++i) {
        ret.push_back(read(i));
    }
    return ret;
}

arma::mat File::read_file(const char *path, int idx) {
    return File::open_file(path).read(idx);
}

void File::write_file(const char *path, Segment segment, const double *data, Format format) {
    auto file = File::create_file(path);
    file.write(std::move(segment), data, format);
}

void File::write(Segment segment, const double *data, Format format) {
    auto ovf_format = format == Format::Text ? OVF_FORMAT_TEXT : OVF_FORMAT_BIN;
    // for some reason, 'ovf_write_segment_8' doesn't take a 'const *', even though it could
    handle_err("write segment", ovf_write_segment_8, segment.ptr.get(), const_cast<double *>(data), ovf_format);
    segments.emplace_back(std::move(segment));
}

}