#include <sstream>
#include <functional>

#include "log_buffer.h"

extern "C" {
#include "cafe.h"
}

using namespace std;

log_buffer::log_buffer(pCafeParam param, size_t buff_sz) : param_(param), buffer_(buff_sz + 1)
{
	char *base = &buffer_.front();
	setp(base, base + buffer_.size() - 1); // -1 to make overflow() easier
}

int log_buffer::sync()
{
	std::ptrdiff_t n = pptr() - pbase();
	ostringstream ost;
	ost.write(pbase(), n);
	pbump(-n);
	cafe_log(param_, ost.str().c_str());
	return 0;
}

log_buffer::int_type log_buffer::overflow(int_type ch)
{
	if (ch != traits_type::eof())
	{
		assert(std::less_equal<char *>()(pptr(), epptr()));
		*pptr() = ch;
		pbump(1);
		sync();
		return ch;
	}

	return traits_type::eof();
}

