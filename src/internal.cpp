#include <gcsa/internal.h>

namespace gcsa
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type DiskIO::block_size;

template<class Element> constexpr size_type ReadBuffer<Element>::READ_BUFFER_SIZE;
template<class Element> constexpr size_type ReadBuffer<Element>::MINIMUM_SIZE;

//------------------------------------------------------------------------------

// Other class variables.

std::atomic<size_type> DiskIO::read_volume(0);
std::atomic<size_type> DiskIO::write_volume(0);

//------------------------------------------------------------------------------

CounterArray::CounterArray() :
  width(8), large_value(sdsl::bits::lo_set[width]),
  total(0)
{
}

CounterArray::CounterArray(size_type n, size_type data_width) :
  data(n, 0, data_width),
  width(data_width), large_value(sdsl::bits::lo_set[width]),
  total(0)
{
}

void
CounterArray::clear()
{
  sdsl::util::clear(this->data);
  sdsl::util::clear(this->large_values);
  this->total = 0;
}

void
CounterArray::swap(CounterArray& another) noexcept
{
  if(this != &another)
  {
    this->data.swap(another.data);
    this->large_values.swap(another.large_values);
    std::swap(this->width, another.width);
    std::swap(this->large_value, another.large_value);
    std::swap(this->total, another.total);
  }
}

//------------------------------------------------------------------------------

} // namespace gcsa
