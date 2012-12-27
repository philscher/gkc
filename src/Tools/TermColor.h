/*
 * =====================================================================================
 *
 *       Filename: TermColor.h
 *
 *    Description: Consol Color support through escape sequences
 *                 (Currently only bash is tested)
 *
 *         Author: Paul P. Hilscher (2012), 
 *
 *        License: MIT
 * =====================================================================================
 */


#ifndef __GKC_TERMCOLOR_H__
#define __GKC_TERMCOLOR_H__

class TermColor

{
  public:
    
    // forground colors
    static const char *cdefault;

    static const char *black  ;
    static const char *red    ;
    static const char *green  ;
    static const char *yellow ;
    static const char *blue   ;
    static const char *magenta;
    static const char *cyan   ;
    static const char *white  ;
    
    // forground light colors
    static const char *lblack  ;
    static const char *lred    ;
    static const char *lgreen  ;
    static const char *lyellow ;
    static const char *lblue   ;
    static const char *lmagenta;
    static const char *lcyan   ;
    static const char *lwhite  ;

    // background colors
    static const char *bcdefault;
    static const char *bblack   ;
    static const char *bred     ;
    static const char *bgreen   ;
    static const char *byellow  ;
    static const char *bblue    ;
    static const char *bmagenta ;
    static const char *bcyan    ;
    static const char *bwhite   ;
    
    // background light colors
    static const char *blblack  ;
    static const char *blred    ;
    static const char *blgreen  ;
    static const char *blyellow ;
    static const char *blblue   ;
    static const char *blmagenta;
    static const char *blcyan   ;
    static const char *blwhite  ;

    // text modifiers
    static const char *reset    ;
    static const char *bright   ;
    static const char *dim      ;
    static const char *underline;
    static const char *blink    ;
    static const char *reverse  ;
    static const char *hidden   ;
};


#endif // __GKC_TERMCOLOR_H__
