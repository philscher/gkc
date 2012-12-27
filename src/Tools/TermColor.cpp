/*
 * =====================================================================================
 *
 *       Filename: TermColor.cpp
 *
 *    Description: Consol Color support through escape sequences
 *                 (Currently only bash is tested)
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: MIT
 * =====================================================================================
 */


#include "TermColor.h"


//// Forground colors
const char *TermColor::cdefault  = "\033[0m";

const char *TermColor::black     = "\033[30m";
const char *TermColor::red       = "\033[31m";
const char *TermColor::green     = "\033[32m";
const char *TermColor::yellow    = "\033[33m";
const char *TermColor::blue      = "\033[34m";
const char *TermColor::magenta   = "\033[35m";
const char *TermColor::cyan      = "\033[36m";
const char *TermColor::white     = "\033[37m";

const char *TermColor::lblack    = "\033[90m";
const char *TermColor::lred      = "\033[91m";
const char *TermColor::lgreen    = "\033[92m";
const char *TermColor::lyellow   = "\033[93m";
const char *TermColor::lblue     = "\033[94m";
const char *TermColor::lmagenta  = "\033[95m";
const char *TermColor::lcyan     = "\033[96m";
const char *TermColor::lwhite    = "\033[97m";


//// Background colors

const char *TermColor::bcdefault = "\033[49m";

const char *TermColor::bblack    = "\033[40m";
const char *TermColor::bred      = "\033[41m";
const char *TermColor::bgreen    = "\033[42m";
const char *TermColor::byellow   = "\033[43m";
const char *TermColor::bblue     = "\033[44m";
const char *TermColor::bmagenta  = "\033[45m";
const char *TermColor::bcyan     = "\033[46m";
const char *TermColor::bwhite    = "\033[47m";

const char *TermColor::blblack   = "\033[100m";
const char *TermColor::blred     = "\033[101m";
const char *TermColor::blgreen   = "\033[102m";
const char *TermColor::blyellow  = "\033[103m";
const char *TermColor::blblue    = "\033[104m";
const char *TermColor::blmagenta = "\033[105m";
const char *TermColor::blcyan    = "\033[106m";
const char *TermColor::blwhite   = "\033[107m";


const char *TermColor::reset     = "\033[0m";
const char *TermColor::bright    = "\033[1m";
const char *TermColor::dim       = "\033[2m";
const char *TermColor::underline = "\033[3m";
const char *TermColor::blink     = "\033[4m";
const char *TermColor::reverse   = "\033[5m";
const char *TermColor::hidden    = "\033[6m";

