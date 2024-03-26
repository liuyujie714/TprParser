#ifndef BYTES_H
#define BYTES_H

#include <stdio.h>
#include <stdexcept>
#include "endianswap.h"
#include <type_traits>
#include <string>

#define TPR_SUCCESS true
#define CPT_SUCCESS true
#define EDR_SUCCESS true
#define TPR_FAILED  false
#define CPT_FAILED  false
#define EDR_FAILED  false
#define MIN(a, b) ((a)<(b)? (a):(b))
#define SAVELEN 512
#define CPT_MAGIC1 171817
#define CPT_MAGIC2 171819

class FileSerializer
{
public:
	FileSerializer(const char *fname, const char *mode) : m_fname(fname)
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

		// note
		fprintf(stderr, "NOTE) Open file %s to %s\n", fname, m_read ? "read" : "write");

		// need endianism swap?
		m_rev = is_litendian();
		msg("data endian= %s\n", m_rev ? "little" : "big");

		// get all buffer
		if (m_read) get_buffer();
	}

	~FileSerializer()
	{
		if (m_buffer) delete [] m_buffer;
		if (fp) fclose(fp);
		fprintf(stderr, "NOTE) End of %s to %s\n", m_fname.c_str(), m_read ? "read" : "write");
	}

	//< get a pointer to file char *buffer
	const char* get_file_buffer(long *fsize) const
	{
		*fsize = m_fsize;
		return m_buffer;
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
				unsigned short tempui = *val;
				if (m_rev) swap2_aligned(&tempui, 1);
				if (fwrite(&tempui, 2, 1, fp) != 1) return TPR_FAILED;
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
			int tempint = *val; // avoid change *val binary order
			if (m_rev) swap4_aligned(&tempint, 1);
			if (fwrite(&tempint, 4, 1, fp) != 1) return TPR_FAILED;
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
			int64_t tempint64 = *val; // avoid change *val
			if (m_rev) swap8_aligned(&tempint64, 1);
			if (fwrite(&tempint64, 8, 1, fp) != 1) return TPR_FAILED;
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
			float tempfloat = *val; // avoid change *val
			if (m_rev) swap4_aligned(&tempfloat, 1);
			if (fwrite(&tempfloat, 4, 1, fp) != 1) return TPR_FAILED;
		}
		return TPR_SUCCESS;
	}

	// read/write double
	bool do_double(double* val) const
	{
		static_assert(sizeof(double) == 8, "sizeof double must be 8");
		if (m_read)
		{
			if (fread(val, 8, 1, fp) != 1) return TPR_FAILED;
			if (m_rev) swap8_aligned(val, 1);
		}
		else
		{
			double tempdouble = *val;
			if (m_rev) swap8_aligned(&tempdouble, 1);
			if (fwrite(&tempdouble, 8, 1, fp) != 1) return TPR_FAILED;
		}
		return TPR_SUCCESS;
	}

	//< read/write float/double according to given prec (4/8)
	bool do_real(void * val, int prec) const
	{
		float	f;
		double	d;
		switch (prec)
		{
			case sizeof(float) :
			{
				if (m_read)
				{
					if (!do_float(&f)) return TPR_FAILED;
					*(static_cast<float*>(val)) = f;
				}
				else
				{
					f = *(static_cast<float*>(val));
					if (!do_float(&f)) return TPR_FAILED;
				}

				break;
			}
			case sizeof(double) :
			{
				if (m_read)
				{
					if (!do_double(&d)) return TPR_FAILED;
					*(static_cast<float*>(val)) = static_cast<float>(d); // double to float
				}
				else
				{
					d = *(static_cast<double*>(val));
					if (!do_double(&d)) return TPR_FAILED;
				}

				break;
			}
			default:
				throw std::runtime_error("Can not support precision= " + std::to_string(prec));
			}
		return TPR_SUCCESS;
	}

	//< read/write us, int, int64_t, float, double, ... in vector with len
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
	bool xdr_string(char* str, int max) const
	{
		int size;
		if (do_int(&size) == TPR_FAILED) return TPR_FAILED;

		// 字符串长度不是4的倍数 {VERSION 2019.6}
		if (size % 4)
		{
			size += 4 - (size % 4); // 满足4的倍数
		}

		size_t ssize = static_cast<size_t>(size); // ignore warning
		if (str && size < max)
		{
			if (fread(str, 1, ssize, fp) != ssize) return TPR_FAILED;
			str[ssize] = '\0';
			return TPR_SUCCESS;
		}
		// size >= max
		else if (str)
		{
			if (fread(str, 1, (size_t)max, fp) != (size_t)max) return TPR_FAILED;
			str[max - 1] = '\0';
			// skip next string
			if (fseek(fp, ssize - max, SEEK_CUR) != 0) return TPR_FAILED;
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

	//< report file pointer position
	long ftell_() const
	{
		return ftell(fp);
	}

	//< fwrite
	size_t fwrite_(const void* buffer, size_t elementsize, size_t count)
	{
		return fwrite(buffer, elementsize, count, fp);
	}

private:
	//< get all binary file buffer
	void get_buffer()
	{
		fseek_(0, SEEK_END); // file end
		long fsize = ftell_();
		m_buffer = new char[fsize];

		fseek_(0, SEEK_SET); // file start
		if (fread(m_buffer, fsize, 1, fp) != 1)
		{
			throw std::runtime_error("Can not read all binary stream to m_buffer");
		}
		// restore
		fseek_(0, SEEK_SET);

		// filesize
		m_fsize = fsize;
	}

private:
	FILE		*fp = nullptr; //< file pointer
	std::string m_fname = {}; //< file name
	bool		m_read = true; //< if read mode
	bool		m_rev = false; //< if Reverse endiannism?
	char		* m_buffer = nullptr; //< all file binary data in char *
	long		m_fsize = 0; // the file size in char
};


#endif // !BYTES_H
