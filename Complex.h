//
//  Simple complex number class
//  D.P. Mitchell 2021/03/02.
//
#include "stdafx.h"
#pragma once
#pragma intrinsic(fabs)

struct Complex {
    double re, im;

            Complex() {}
            Complex(double fr) : re(fr), im(0.0) {}
            Complex(double fr, double fi) : re(fr), im(fi) {}
    int     operator==(Complex c) { return re == c.re && im == c.im; }
    Complex operator+(Complex c) const { return Complex(re + c.re, im + c.im); }
    Complex operator-(Complex c) const { return Complex(re - c.re, im - c.im); }
    Complex operator*(Complex c) const { return Complex(re*c.re - im*c.im, re*c.im + im*c.re); }
    Complex operator*(double f) const  { return Complex(re*f, im*f); }
    Complex operator/(double f) const  { return Complex(re/f, im/f); }
    Complex operator/(Complex c) {
        double fDenom, Q;

        if (fabs(c.re) < fabs(c.im)) {              // overflow-avoiding method
            Q = c.re/c.im;
            fDenom = c.re*Q + c.im;
            return Complex((re*Q + im)/fDenom, (im*Q - re)/fDenom);
        } else {
            Q = c.im/c.re;
            fDenom = c.im*Q + c.re;
            return Complex((im*Q + re)/fDenom, (im - re*Q)/fDenom);
        }
    }
    double  AbsSqr() const { return re*re + im*im; }
};