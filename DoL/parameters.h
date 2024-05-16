#ifndef DOL_PARAMETERS_H_INCLUDED
#define DOL_PARAMETERS_H_INCLUDED

#include <string>
#include "cmd_line.h"
#include <fstream>

struct param_t
{
  int Tmax = 10000;
  double forget_rate = 0.9;
  double learning_rate = 0.1; 
  double updaterate = 0.2;
  double mort = 0.006; //0.006
  double birthrate = 0.08; //0.2 //6, 0.2, 0.06
  double f1 = 0.5;
  int seed = 1;
  std::string outdir;



};

cmd::cmd_line_parser config_file_parser(const std::string& config)
{
  std::ifstream is(config);
  if (!is) throw cmd::parse_error("can't open config file");
  std::string str((std::istreambuf_iterator<char>(is)), std::istreambuf_iterator<char>());
  return cmd::cmd_line_parser(str);
}



#endif