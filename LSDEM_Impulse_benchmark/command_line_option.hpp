#ifndef COMMAND_LINE_OPTION_HPP_
#define COMMAND_LINE_OPTION_HPP_
int find_arg_idx(int argc, char** argv, const char* option) {
    for(int i = 1; i < argc; ++i) {
        if(strcmp(argv[i], option) == 0) return i;
    }
    return -1;
}

int find_int_arg(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);
    if(iplace >= 0 && iplace < argc - 1) return std::stoi(argv[iplace + 1]);
    return default_value;
}

double find_double_arg(int argc, char** argv, const char* option, double default_value) {
    int iplace = find_arg_idx(argc, argv, option);
    if(iplace >= 0 && iplace < argc - 1) return std::stod(argv[iplace + 1]);
    return default_value;
}

char* find_string_option(int argc, char** argv, const char* option, char* default_value) {
    int iplace = find_arg_idx(argc, argv, option);
    if(iplace >= 0 && iplace < argc - 1) return argv[iplace + 1];
    return default_value;
}

int find_string_option(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);
    if(iplace >= 0 && iplace < argc - 1){
      if(strcmp(argv[iplace + 1], "flexible") == 0) { return 0; }
      else if(strcmp(argv[iplace + 1], "rigid") == 0) { return 1; }
      else if(strcmp(argv[iplace + 1], "cylinder") == 0) { return 2; }
      else if(strcmp(argv[iplace + 1], "lateral") == 0) { return 3; }
    }
    return default_value;
}
#endif /*COMMAND_LINE_OPTION_HPP_*/
