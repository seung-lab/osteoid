#ifndef _OSTD_LIB_HXX_
#define _OSTD_LIB_HXX_

namespace fastosteoid::lib {

template <typename T>
T ctoi(const unsigned char* buf, const uint64_t idx = 0);

template <>
int64_t ctoi(const unsigned char* buf, const uint64_t idx) {
	int64_t x = 0;
	x |= static_cast<uint64_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint64_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint64_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint64_t>(buf[idx + 3]) << 24;
	x |= static_cast<uint64_t>(buf[idx + 4]) << 32;
	x |= static_cast<uint64_t>(buf[idx + 5]) << 40;
	x |= static_cast<uint64_t>(buf[idx + 6]) << 48;
	x |= static_cast<uint64_t>(buf[idx + 7]) << 56;
	return x;
}

template <>
uint64_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint64_t x = 0;
	x |= static_cast<uint64_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint64_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint64_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint64_t>(buf[idx + 3]) << 24;
	x |= static_cast<uint64_t>(buf[idx + 4]) << 32;
	x |= static_cast<uint64_t>(buf[idx + 5]) << 40;
	x |= static_cast<uint64_t>(buf[idx + 6]) << 48;
	x |= static_cast<uint64_t>(buf[idx + 7]) << 56;
	return x;
}

template <>
int32_t ctoi(const unsigned char* buf, const uint64_t idx) {
	int32_t x = 0;
	x |= static_cast<uint32_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint32_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint32_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint32_t>(buf[idx + 3]) << 24;
	return x;
}

template <>
uint32_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint32_t x = 0;
	x |= static_cast<uint32_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint32_t>(buf[idx + 1]) << 8;
	x |= static_cast<uint32_t>(buf[idx + 2]) << 16;
	x |= static_cast<uint32_t>(buf[idx + 3]) << 24;
	return x;
}

template <>
int16_t ctoi(const unsigned char* buf, const uint64_t idx) {
	int16_t x = 0;
	x |= static_cast<uint16_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint16_t>(buf[idx + 1]) << 8;
	return x;
}

template <>
uint16_t ctoi(const unsigned char* buf, const uint64_t idx) {
	uint16_t x = 0;
	x |= static_cast<uint16_t>(buf[idx + 0]) << 0;
	x |= static_cast<uint16_t>(buf[idx + 1]) << 8;
	return x;
}

template <>
uint8_t ctoi(const unsigned char* buf, const uint64_t idx) {
	return static_cast<uint8_t>(buf[idx]);
}

template <>
int8_t ctoi(const unsigned char* buf, const uint64_t idx) {
	return static_cast<int8_t>(buf[idx]);
}

template <typename T>
T ctoi(const std::span<const uint8_t>& buf, const uint64_t idx) {
	return ctoi<T>(buf.data(), idx);
}

uint64_t ctoid(
	const unsigned char* buf, const uint64_t idx, const int byte_width
) {
	uint64_t val = 0;
	for (int i = 0; i < byte_width; i++) {
		val |= (buf[idx + i] << (i*8));
	}
	return val;
}

uint64_t ctoid(
	const std::vector<unsigned char>& buf,
	const uint64_t idx, const int byte_width
) {
	return ctoid(buf.data(), idx, byte_width);
}

uint64_t ctoid(
	const std::span<const unsigned char>& buf,
	const uint64_t idx, const int byte_width
) {
	return ctoid(buf.data(), idx, byte_width);
}

};

#endif