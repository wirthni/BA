#include <cstdlib>
#include <cxxabi.h>   // demangle
#include <typeinfo>
#include <array>

#include <iostream>   // cout used only in the test

using namespace std;

std::string demangle(const char* symbol) {
  int     status;

  char *realname;
  realname = abi::__cxa_demangle(symbol, NULL, NULL, &status);

  std::string demangled = (status==0)
    ? realname : symbol;
  std::free(realname);

  return demangled;
}

int main(void){ 

  // typeid
  std::array<int,17>  u;
  const std::type_info  &ti = typeid(u);

  std::cout << ti.name() << "\t=> " << demangle(ti.name()) << '\n';

}
