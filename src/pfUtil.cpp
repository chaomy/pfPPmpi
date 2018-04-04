/*
 * @Author: chaomy
 * @Date:   2017-11-04 14:53:19
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-04-04 13:33:37
 */

#include "pfHome.h"

void pfUtil::split(const string& s, const char* delim, vector<string>& v) {
  // first duplicate the original string and return a char pointer then free the
  // memory
  char* dup = strdup(s.c_str());
  char* token = strtok(dup, delim);
  while (token != NULL) {
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
}

void pfHome::outMkdir(string mdir){
    struct stat buf;
    if ((stat(mdir.c_str(), &buf) == 0) == 1)
        printf("kmcOutput already exists.\n");
    else
        mkdir(mdir.c_str(), S_IRWXU); 
}