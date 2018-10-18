// -*- C++ -*-
//
//--------------------------------------------------------------------------------
//
//                              Michael A.G. Aivazis
//                       California Institute of Technology
//                       (C) 1998-2003  All Rights Reserved
//
// <LicenseText>
//
//--------------------------------------------------------------------------------
//

#if !defined(journal_manip_decls_h)
#define journal_manip_decls_h

namespace journal {

    // manipulators with no arguments
    inline Diagnostic & operator<< (Diagnostic &, Diagnostic & (*)(Diagnostic &));

    // manipulators with one argument
    template <typename arg1_t> class manip_1;
    template <typename arg1_t>
    Diagnostic & operator << (Diagnostic &, manip_1<arg1_t>);

    // manipulators with two arguments
    template <typename arg1_t, typename arg2_t> class manip_2;
    template <typename arg1_t, typename arg2_t>
    Diagnostic & operator << (Diagnostic &, manip_2<arg1_t, arg2_t>);

    // manipulators with three arguments
    template <typename arg1_t, typename arg2_t, typename arg3_t>
    class manip_3;
    template <typename arg1_t, typename arg2_t, typename arg3_t>
    Diagnostic & operator << (Diagnostic &, manip_3<arg1_t, arg2_t, arg3_t>);

}


template <typename arg1_t>
class journal::manip_1 {

    friend Diagnostic & operator<< <> (Diagnostic &, manip_1<arg1_t>);

// types
public:
    typedef Diagnostic & (*factory_t)(Diagnostic &, arg1_t);

// meta-methods
public:

    manip_1(factory_t f, arg1_t arg1);

// data
private:

    factory_t _f;
    arg1_t _arg1;
};

template <typename arg1_t, typename arg2_t>
class journal::manip_2 {

    friend Diagnostic & operator<< <> (Diagnostic &, manip_2<arg1_t, arg2_t>);

// types
public:
    typedef Diagnostic & (*factory_t)(Diagnostic &, arg1_t, arg2_t);

// meta-methods
public:

    manip_2(factory_t f, arg1_t arg1, arg2_t arg2);

// data
private:

    factory_t _f;
    arg1_t _arg1;
    arg2_t _arg2;
};

template <typename arg1_t, typename arg2_t, typename arg3_t>
class journal::manip_3 {

    friend Diagnostic & operator<< <> (Diagnostic &, manip_3<arg1_t, arg2_t, arg3_t>);

// types
public:
    typedef Diagnostic & (*factory_t)(Diagnostic &, arg1_t, arg2_t, arg3_t);

// meta-methods
public:

    manip_3(factory_t f, arg1_t arg1, arg2_t arg2, arg3_t arg3);

// data
private:

    factory_t _f;
    arg1_t _arg1;
    arg2_t _arg2;
    arg3_t _arg3;
};

#include "manip-decls.icc"

#endif

// version
// $Id$

// End of file
