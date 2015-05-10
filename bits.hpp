#pragma once
#include <cstdint>

//! 16bitの数値ビットを1つ飛びに変換する
uint32_t SeparateBits1(uint16_t val);
//! 32bitの1つ飛び16bit数値を元に戻す
uint32_t ComposeBits1(uint32_t x);
//! 10bitまでの数値を2つ飛びに変換する
uint32_t SeparateBits2(uint32_t x);
//! 32bitの2つ飛び10bit数値を元に戻す
uint32_t ComposeBits2(uint32_t x);

