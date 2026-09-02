#include <gcsa/internal.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>

using namespace gcsa;

static void require(bool value)
{
  if(!value) { std::abort(); }
}

int main()
{
  const std::string name = "test_internal_buffers_" +
    std::to_string(static_cast<unsigned long long>(getpid())) + ".bin";

  {
    WriteBuffer<size_type> output(name, 33);
    require(output.buffer_size == 33 / sizeof(size_type));
    require(output.buffer.capacity() == output.buffer_size);
    for(size_type i = 0; i < 17; i++) { output.push_back(i * i); }
  }

  {
    ReadBuffer<size_type> input;
    input.open(name, 33);
    require(input.buffer_elements == 33 / sizeof(size_type));
    require(input.minimum_elements == 2);
    require(input.size() == 17);
    for(size_type i = 0; i < input.size(); i++) { require(input[i] == i * i); }
    input.close();
  }

  {
    const std::string bytes_name = name + ".bytes";
    WriteBuffer<uint8_t> output(bytes_name, 33);
    require(output.buffer_size == 33);
    output.push_back(1);
    output.close();
    std::remove(bytes_name.c_str());
  }

  std::remove(name.c_str());
  return 0;
}
