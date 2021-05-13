// STATUS BAR  
// Modified
// http://www.rosshemsley.co.uk/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/

#include <iomanip>

static inline void loadBar(unsigned int x, unsigned int n, unsigned int w = 50)
{

  //Make sure that load bar is displayed on next line and don’t delete current console line during first run
  static bool firstCall = true;
  if (firstCall)
  {
    std::cout << std::endl;
    firstCall = false;
  }

  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  float ratio_ = x / (float)n;
  int c = int(ratio_ * w);

  std::cout << std::setw(3) << (int)(ratio_ * 100) << "% [";
  for (int x = 0; x < c; x++)
  {
    std::cout << "\033[1;32m|";
  }
  for (unsigned int x = c; x < w; x++)
  {
    std::cout << " ";
  }
  std::cout << "\033[0m]\r" << std::flush;
}
