#include "bits.hpp"
#include <cassert>

uint32_t SeparateBits1(uint16_t val) {
	uint32_t x = val;
	x = (x | (x<<8)) & 0x00ff00ff;
	x = (x | (x<<4)) & 0x0f0f0f0f;
	x = (x | (x<<2)) & 0x33333333;
	x = (x | (x<<1)) & 0x55555555;
	return x;
}
uint32_t ComposeBits1(uint32_t x) {
	x = (x | (x>>1)) & 0x33333333;
	x = (x | (x>>2)) & 0x0f0f0f0f;
	x = (x | (x>>4)) & 0x00ff00ff;
	x = (x | (x>>8)) & 0x0000ffff;
	return x;
}
uint32_t SeparateBits2(uint32_t x) {
	assert(x < (1<<10));				// 00000000 : 00000000 : 00000012 : 3456789A
	x = (x | (x<<16)) & 0x030000ff;		// 00000012 : 00000000 : 00000000 : 3456789A
	x = (x | (x<<8)) & 0x0300f00f;		// 00000012 : 00000000 : 34560000 : 0000789A
	x = (x | (x<<4)) & 0x030c30c3;		// 00000012 : 00003400 : 00560000 : 7800009A
	x = (x | (x<<2)) & 0x09249249;		// 00001002 : 00300400 : 50060070 : 0800900A
	return x;
}
uint32_t ComposeBits2(uint32_t x) {
	x = (x | (x>>2)) & 0x030c30c3;
	x = (x | (x>>4)) & 0x0300f00f;
	x = (x | (x>>8)) & 0x030000ff;
	x = (x | (x>>16)) & 0x000003ff;
	return x;
}

