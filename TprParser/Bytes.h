#ifndef BYTES_H
#define BYTES_H

#include <stdio.h>
#include <stdexcept>
#include "endianswap.h"
#include <type_traits>
#include <string>

#define TPR_SUCCESS true
#define TPR_FAILED  false
#define MIN(a, b) ((a)<(b)? (a):(b))
#define SAVELEN 512

class TpxSerializer
{
public:
	TpxSerializer(const char *fname, const char *mode)
	{
		fp = fopen(fname, mode);
		if (!fp)
		{
			throw std::runtime_error("Can not open/write file: " + std::string(fname));
		}
		// is read model
		switch (mode[0])
		{
		case 'r':
			m_read = true;
			break;
		case 'w':
			m_read = false;
			break;
		}

		// need endianism swap?
		m_rev = is_litendian();
		if (m_rev)
		{
			msg("is_litendian\n");
		}
	}

	~TpxSerializer()
	{
		if (fp) fclose(fp);
	}

	// if is little endian
	bool is_litendian() const
	{
		union 
		{
			int a;
			char b;
		} u;
		u.a = 1;
		return u.b == 1;
	}

	// read/write bool, return TPR_SUCCESS if succeed
	bool do_bool(bool* val, int vergen = 26) const
	{
		if (m_read)
		{
			// read 1 byte
			if (vergen >= 27)
			{
				if (fread(val, 1, 1, fp) != 1) return TPR_FAILED;
			}
			else
			{
				int tempint = 0;
				if (fread(&tempint, 4, 1, fp) != 1) return TPR_FAILED;
				if (m_rev) swap4_aligned(&tempint, 1);
				*val = (tempint != 0); // return bool
			}
		}
		else
		{
			// write 1 byte
			if (vergen >= 27)
			{
				if (fwrite(val, 1, 1, fp) != 1) return TPR_FAILED;
			}
			else
			{
				// bool to int
				int tempint = static_cast<int>(*val);
				if (m_rev) swap4_aligned(&tempint, 1);
				if (fwrite(&tempint, 4, 1, fp) != 1) return TPR_FAILED;
			}
		}
		return TPR_SUCCESS;
	}

	// read/write unsigned short, return TPR_SUCCESS if succeed
	// actually read int if vergen < 27 and convert to unsigned short
	bool do_ushort(unsigned short * val, int vergen = 26) const
	{
		static_assert(sizeof(unsigned short) == 2, "sizeof unsigned short must be 2");
		if (m_read)
		{
			// for gmx2020
			if (vergen >= 27)
			{
				if (fread(val, 2, 1, fp) != 1) return TPR_FAILED;
				if (m_rev) swap2_aligned(val, 1);
			}
			else
			{
				int temp;
				if (fread(&temp, 4, 1, fp) != 1) return TPR_FAILED;
				if (m_rev) swap4_aligned(&temp, 1);
				*val = static_cast<unsigned short>(temp);
			}
		}
		else
		{
			// for gmx2020
			if (vergen >= 27)
			{
				if (m_rev) swap2_aligned(val, 1);
				if (fwrite(val, 2, 1, fp) != 1) return TPR_FAILED;
			}
			else
			{
				int temp = static_cast<int>(*val);
				if (m_rev) swap4_aligned(&temp, 1);
				if (fwrite(&temp, 4, 1, fp) != 1) return TPR_FAILED;
			}
		}

		return TPR_SUCCESS;
	}

	// read/write unsigned char, return TPR_SUCCESS if succeed
	// actually read unsigned int if vergen < 27 and convert to unsigned short
	bool do_uchar(unsigned char * val, int vergen = 26) const
	{
		static_assert(sizeof(unsigned char) == 1, "sizeof unsigned char must be 1");
		if (m_read)
		{
			// for gmx>=2020, only read 1 byte
			if (vergen >= 27)
			{
				if (fread(val, 1, 1, fp) != 1) return TPR_FAILED;
			}
			else
			{
				int temp;
				if (fread(&temp, 4, 1, fp) != 1) return TPR_FAILED;
				if (m_rev) swap4_aligned(&temp, 1);
				*val = static_cast<unsigned char>(temp);
			}
		}
		else
		{
			// for gmx>=2020, only write 1 byte
			if (vergen >= 27)
			{
				if (fwrite(val, 1, 1, fp) != 1) return TPR_FAILED;
			}
			else
			{
				int temp = static_cast<int>(*val);
				if (m_rev) swap4_aligned(&temp, 1);
				if (fwrite(&temp, 4, 1, fp) != 1) return TPR_FAILED;
			}
		}

		return TPR_SUCCESS;
	}

	// read/write int32, return TPR_SUCCESS if succeed
	bool do_int(int *val) const
	{
		static_assert(sizeof(int) == 4, "sizeof int must be 4");
		if (m_read)
		{
			if (fread(val, 4, 1, fp) != 1) return TPR_FAILED;
			if (m_rev) swap4_aligned(val, 1);
		}
		else
		{
			if (m_rev) swap4_aligned(val, 1);
			if (fwrite(val, 4, 1, fp) != 1) return TPR_FAILED;
		}
		return TPR_SUCCESS;
	}

	// read/write int64
	bool do_int64(int64_t * val) const
	{
		static_assert(sizeof(int64_t) == 8, "sizeof int64_t must be 8");
		if (m_read)
		{
			if (fread(val, 8, 1, fp) != 1) return TPR_FAILED;
			if (m_rev) swap8_aligned(val, 1);
		}
		else
		{
			if (m_rev) swap8_aligned(val, 1);
			if (fwrite(val, 8, 1, fp) != 1) return TPR_FAILED;
		}
		return TPR_SUCCESS;
	}

	// read/write float
	bool do_float(float *val) const
	{
		static_assert(sizeof(float) == 4, "sizeof float must be 4");
		if (m_read)
		{
			if (fread(val, 4, 1, fp) != 1) return TPR_FAILED;
			if (m_rev) swap4_aligned(val, 1);
		}
		else
		{
			if (m_rev) swap4_aligned(val, 1);
			if (fwrite(val, 4, 1, fp) != 1) return TPR_FAILED;
		}
		return TPR_SUCCESS;
	}

	// read/write double
	bool do_double(double* val) const
	{
		static_assert(sizeof(double) == 8, "sizeof double must be 4");
		if (m_read)
		{
			if (fread(val, 8, 1, fp) != 1) return TPR_FAILED;
			if (m_rev) swap8_aligned(val, 1);
		}
		else
		{
			if (m_rev) swap8_aligned(val, 1);
			if (fwrite(val, 8, 1, fp) != 1) return TPR_FAILED;
		}
		return TPR_SUCCESS;
	}

	//< only read float/double according to given prec (4/8)
	bool do_real(void * val, int prec) const
	{
		float	f;
		double	d;
		switch (prec)
		{
			case sizeof(float) :
			{
				if (!do_float(&f)) return TPR_FAILED;
				*(static_cast<float*>(val)) = f;
				break;
			}
			case sizeof(double) :
			{
				if (!do_double(&d)) return TPR_FAILED;
				*(static_cast<float*>(val)) = static_cast<float>(d); // double to float
				break;
			}
			default:
				throw std::runtime_error("Can not support precision= " + std::to_string(prec));
			}
		return TPR_SUCCESS;
	}

	//< read float in len vector
	template<typename T>
	bool do_vector(T* arr, int len, int prec = 4, int vergen = 26) const
	{
		for (int i = 0; i < len; i++)
		{
			if constexpr (std::is_same_v<T, unsigned char>)
			{
				if (!do_uchar(&arr[i], vergen)) return TPR_FAILED;
			}
			else if constexpr (std::is_same_v<T, int>)
			{
				if (!do_int(&arr[i])) return TPR_FAILED;
			}
			else if constexpr (std::is_same_v <T, int64_t>)
			{
				if (!do_int64(&arr[i])) return TPR_FAILED;
			}
			else if constexpr (std::is_same_v <T, float>)
			{
				if (!do_real(&arr[i], prec)) return TPR_FAILED;
			}
			else if constexpr (std::is_same_v <T, double>)
			{
				if (!do_double(&arr[i])) return TPR_FAILED;
			}
			else
			{
				throw std::runtime_error("Unsupport type for do_vector\n");
			}
		}
		return TPR_SUCCESS;
	}

	// Reads in a string by first reading an integer containing the
	// string's length, then reading in the string itself and storing
	// it in str. If the length is greater than max, it is truncated
	// and the rest of the string is skipped in the file
	bool tpr_string(char* str, int max) const
	{
		int size;
		if (do_int(&size) == TPR_FAILED) return TPR_FAILED;

		// 字符串长度不是4的倍数 {VERSION 2019.6}
		if (size % 4)
		{
			size += 4 - (size % 4); // 满足4的倍数
		}
		if (str && size <= max)
		{
			if (fread(str, 1, size, fp) != size) return TPR_FAILED;
			str[size] = '\0';
			return TPR_SUCCESS;
		}
		// size > max
		else if (str)
		{
			if (fread(str, 1, max, fp) != max) return TPR_FAILED;
			// skip next string
			if (fseek(fp, size - max, SEEK_CUR) != 0) return TPR_FAILED;
			str[max] = '\0';
			return TPR_SUCCESS;
		}
		else
		{
			// skip all string and don not store
			if (fseek(fp, size, SEEK_CUR) != 0) return TPR_FAILED;
			return TPR_SUCCESS;
		}
	}

	// save string to saveloc
	bool tpr_save_string(char* saveloc, int genversion) 
	{
		int			i;
		char		buf[MAX_LEN];

		// for gmx>=2020
		if (genversion >= 27)
		{
			int64_t len;
			if (!do_int64(&len)) return TPR_FAILED;
			if (fread(buf, 1, (size_t)(len), fp) != (size_t)(len)) return TPR_FAILED;
			for (i = 0; i < MIN(len, (SAVELEN - 1)); i++) {
				saveloc[i] = buf[i];
			}
			saveloc[i] = '\0';
		}
		else
		{
			int len;
			// first byte not used
			if (!do_int(&len)) return TPR_FAILED;
			// actually len
			if (!do_int(&len)) return TPR_FAILED;
			if (len % 4) len += 4 - len % 4; // 字节对齐
			if (fread(buf, 1, len, fp) != len) return TPR_FAILED;
			for (i = 0; i < MIN(len, (SAVELEN - 1)); i++) {
				saveloc[i] = buf[i];
			}
			saveloc[i] = '\0';
		}

		return TPR_SUCCESS;
	}

	int fseek_(long offset, int orig)
	{
		return fseek(fp, offset, orig);
	}

private:
	FILE		*fp = nullptr; //< file pointer
	bool		m_read = true; //< if read mode
	bool		m_rev = false; //< if Reverse endiannism?
};


#endif // !BYTES_H
