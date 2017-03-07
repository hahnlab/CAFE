#ifndef LOG_BUFFER_H_BFB0BCB2_DD7C_421B_8BCC_039B98DC0072
#define LOG_BUFFER_H_BFB0BCB2_DD7C_421B_8BCC_039B98DC0072

#include <streambuf>
#include <vector>

extern "C" {
#include <family.h>
}

class log_buffer : public std::streambuf
{
public:
	explicit log_buffer(pCafeParam param, std::size_t buff_sz = 256);
	virtual ~log_buffer()
	{
		sync();
	}
private:
	int_type overflow(int_type ch);
	int sync();

	// copy ctor and assignment not implemented;
	// copying not allowed
	log_buffer(const log_buffer &);
	log_buffer &operator= (const log_buffer &) { return *this; }

private:
	pCafeParam param_;
	std::vector<char> buffer_;
};

#endif
