#pragma once
#ifndef NV_CORE_RADIXSORT_H
#define NV_CORE_RADIXSORT_H

// Based on Pierre Terdiman's and Michael Herf's source code.
// http://www.codercorner.com/RadixSortRevisited.htm
// http://www.stereopsis.com/radix.html

#include <vector>
#include "nvcore.h"
#include "Array.h"

namespace nv
{

class RadixSort
{
public:
	// Constructor/Destructor
	RadixSort();
	~RadixSort();

	RadixSort &sort(const float *input, uint count);
	inline RadixSort &sort(const std::vector<float> &input)
	{
		return sort(input.data(), input.size());
	}

	// Access to results. m_ranks is a list of indices in sorted order, i.e. in the order you may further process your data
	inline const uint *ranks() const
	{
		nvDebugCheck(m_validRanks);
		return m_ranks;
	}
	inline uint *ranks()
	{
		nvDebugCheck(m_validRanks);
		return m_ranks;
	}

private:
	uint m_size;
	uint *m_ranks;
	uint *m_ranks2;
	bool m_validRanks;

	// Internal methods
	template <typename T> void insertionSort(const T *input, uint count);
	template <typename T> void radixSort(const T *input, uint count);
};

} // nv namespace



#endif // NV_CORE_RADIXSORT_H
