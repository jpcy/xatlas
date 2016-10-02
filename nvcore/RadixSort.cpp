// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>

#include "RadixSort.h"

#include "Utils.h"

#include <string.h> // memset

using namespace nv;

static inline void FloatFlip(uint32 &f)
{
	//uint32 mask = -int32(f >> 31) | 0x80000000; // Michael Herf.
	int32 mask = (int32(f) >> 31) | 0x80000000; // Warren Hunt, Manchor Ko.
	f ^= mask;
}

static inline void IFloatFlip(uint32 &f)
{
	uint32 mask = ((f >> 31) - 1) | 0x80000000; // Michael Herf.
	//uint32 mask = (int32(f ^ 0x80000000) >> 31) | 0x80000000; // Warren Hunt, Manchor Ko. @@ Correct, but fails in release on gcc-4.2.1
	f ^= mask;
}

template<typename T>
void createHistograms(const T *buffer, uint count, uint *histogram)
{
	const uint bucketCount = sizeof(T); // (8 * sizeof(T)) / log2(radix)
	// Init bucket pointers.
	uint *h[bucketCount];
	for (uint i = 0; i < bucketCount; i++) {
		h[i] = histogram + 256 * i;
	}
	// Clear histograms.
	memset(histogram, 0, 256 * bucketCount * sizeof(uint));
	// @@ Add support for signed integers.
	// Build histograms.
	const uint8 *p = (const uint8 *)buffer;  // @@ Does this break aliasing rules?
	const uint8 *pe = p + count * sizeof(T);
	while (p != pe) {
		h[0][*p++]++, h[1][*p++]++, h[2][*p++]++, h[3][*p++]++;
		if (bucketCount == 8) h[4][*p++]++, h[5][*p++]++, h[6][*p++]++, h[7][*p++]++;
	}
}

RadixSort::RadixSort() : m_size(0), m_ranks(NULL), m_ranks2(NULL), m_validRanks(false)
{
}

RadixSort::~RadixSort()
{
	// Release everything
	free(m_ranks2);
	free(m_ranks);
}

template <typename T> inline void RadixSort::insertionSort(const T *input, uint count)
{
	if (!m_validRanks) {
		m_ranks[0] = 0;
		for (uint i = 1; i != count; ++i) {
			int rank = m_ranks[i] = i;
			uint j = i;
			while (j != 0 && input[rank] < input[m_ranks[j - 1]]) {
				m_ranks[j] = m_ranks[j - 1];
				--j;
			}
			if (i != j) {
				m_ranks[j] = rank;
			}
		}
		m_validRanks = true;
	} else {
		for (uint i = 1; i != count; ++i) {
			int rank = m_ranks[i];
			uint j = i;
			while (j != 0 && input[rank] < input[m_ranks[j - 1]]) {
				m_ranks[j] = m_ranks[j - 1];
				--j;
			}
			if (i != j) {
				m_ranks[j] = rank;
			}
		}
	}
}

template <typename T> inline void RadixSort::radixSort(const T *input, uint count)
{
	const uint P = sizeof(T); // pass count
	// Allocate histograms & offsets on the stack
	uint histogram[256 * P];
	uint *link[256];
	createHistograms(input, count, histogram);
	// Radix sort, j is the pass number (0=LSB, P=MSB)
	for (uint j = 0; j < P; j++) {
		// Pointer to this bucket.
		const uint *h = &histogram[j * 256];
		const uint8 *inputBytes = (const uint8 *)input; // @@ Is this aliasing legal?
		inputBytes += j;
		if (h[inputBytes[0]] == count) {
			// Skip this pass, all values are the same.
			continue;
		}
		// Create offsets
		link[0] = m_ranks2;
		for (uint i = 1; i < 256; i++) link[i] = link[i - 1] + h[i - 1];
		// Perform Radix Sort
		if (!m_validRanks) {
			for (uint i = 0; i < count; i++) {
				*link[inputBytes[i * P]]++ = i;
			}
			m_validRanks = true;
		} else {
			for (uint i = 0; i < count; i++) {
				const uint idx = m_ranks[i];
				*link[inputBytes[idx * P]]++ = idx;
			}
		}
		// Swap pointers for next pass. Valid indices - the most recent ones - are in m_ranks after the swap.
		swap(m_ranks, m_ranks2);
	}
	// All values were equal, generate linear ranks.
	if (!m_validRanks) {
		for (uint i = 0; i < count; i++) {
			m_ranks[i] = i;
		}
		m_validRanks = true;
	}
}

RadixSort &RadixSort::sort(const float *input, uint count)
{
	if (input == NULL || count == 0) return *this;
	// Resize lists if needed
	if (count != m_size) {
		if (count > m_size) {
			m_ranks2 = (uint *)realloc(m_ranks2, sizeof(uint) * count);
			m_ranks = (uint *)realloc(m_ranks, sizeof(uint) * count);
		}
		m_size = count;
		m_validRanks = false;
	}
	if (count < 32) {
		insertionSort(input, count);
	} else {
		// @@ Avoid touching the input multiple times.
		for (uint i = 0; i < count; i++) {
			FloatFlip((uint32 &)input[i]);
		}
		radixSort<uint32>((const uint32 *)input, count);
		for (uint i = 0; i < count; i++) {
			IFloatFlip((uint32 &)input[i]);
		}
	}
	return *this;
}
