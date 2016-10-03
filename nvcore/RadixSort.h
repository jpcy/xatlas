#pragma once
#ifndef NV_CORE_RADIXSORT_H
#define NV_CORE_RADIXSORT_H

// Based on Pierre Terdiman's and Michael Herf's source code.
// http://www.codercorner.com/RadixSortRevisited.htm
// http://www.stereopsis.com/radix.html

#include <vector>
#include "nvcore.h"
#include "Debug.h"

namespace nv
{

class RadixSort
{
public:
	// Constructor/Destructor
	RadixSort();
	~RadixSort();

	RadixSort &sort(const float *input, uint32_t count);
	inline RadixSort &sort(const std::vector<float> &input)
	{
		return sort(input.data(), input.size());
	}

	// Access to results. m_ranks is a list of indices in sorted order, i.e. in the order you may further process your data
	inline const uint32_t *ranks() const
	{
		nvDebugCheck(m_validRanks);
		return m_ranks;
	}
	inline uint32_t *ranks()
	{
		nvDebugCheck(m_validRanks);
		return m_ranks;
	}

private:
	uint32_t m_size;
	uint32_t *m_ranks;
	uint32_t *m_ranks2;
	bool m_validRanks;

	// Internal methods
	template <typename T> void insertionSort(const T *input, uint32_t count);
	template <typename T> void radixSort(const T *input, uint32_t count);
};

} // nv namespace



#endif // NV_CORE_RADIXSORT_H
