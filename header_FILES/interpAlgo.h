#pragma once

#include"integerMath.h"
#include"helperF.h"
#include"int128g.hpp"

int newtonInterp(LONG *x,LONG *y,const int n,const LONG p,recint P);
int newtonInterp2(LONG* x,
    LONG* y,
    const int n,
    const LONG p);
