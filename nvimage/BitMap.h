// This code is in the public domain -- castanyo@yahoo.es

#pragma once
#ifndef NV_IMAGE_BITMAP_H
#define NV_IMAGE_BITMAP_H

#include "nvcore/nvcore.h"

namespace nv
{
/// Bit map. This should probably be called BitImage.
class BitMap
{
public:
	BitMap() : m_width(0), m_height(0) {}
	BitMap(uint32_t w, uint32_t h) : m_width(w), m_height(h), m_bitArray(w * h) {}

	uint32_t width() const
	{
		return m_width;
	}
	uint32_t height() const
	{
		return m_height;
	}

	void resize(uint32_t w, uint32_t h, bool initValue);

	bool bitAt(uint32_t x, uint32_t y) const
	{
		nvDebugCheck(x < m_width && y < m_height);
		return m_bitArray.bitAt(y * m_width + x);
	}

	void setBitAt(uint32_t x, uint32_t y)
	{
		nvDebugCheck(x < m_width && y < m_height);
		m_bitArray.setBitAt(y * m_width + x);
	}

	void clearAll()
	{
		m_bitArray.clearAll();
	}

private:

	uint32_t m_width;
	uint32_t m_height;
	BitArray m_bitArray;

};

} // nv namespace

#endif // NV_IMAGE_BITMAP_H
