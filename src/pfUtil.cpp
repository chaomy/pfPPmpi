/*
 * @Author: chaomy
 * @Date:   2017-11-04 14:53:19
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-06-26 16:20:02
 */

#include "pfHome.h"
#include "pfIO.h"

void pfHome::pfIO::outMkdir(string mdir) {
  struct stat buf;
  if ((stat(mdir.c_str(), &buf) == 0) == 1)
    printf("kmcOutput already exists.\n");
  else
    mkdir(mdir.c_str(), S_IRWXU);
}