// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME SliceStructDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "slice_struct.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new__pfp(void *p = nullptr);
   static void *newArray__pfp(Long_t size, void *p);
   static void delete__pfp(void *p);
   static void deleteArray__pfp(void *p);
   static void destruct__pfp(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::_pfp*)
   {
      ::_pfp *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::_pfp >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("_pfp", ::_pfp::Class_Version(), "slice_struct.h", 7,
                  typeid(::_pfp), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::_pfp::Dictionary, isa_proxy, 4,
                  sizeof(::_pfp) );
      instance.SetNew(&new__pfp);
      instance.SetNewArray(&newArray__pfp);
      instance.SetDelete(&delete__pfp);
      instance.SetDeleteArray(&deleteArray__pfp);
      instance.SetDestructor(&destruct__pfp);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::_pfp*)
   {
      return GenerateInitInstanceLocal(static_cast<::_pfp*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::_pfp*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new__slice(void *p = nullptr);
   static void *newArray__slice(Long_t size, void *p);
   static void delete__slice(void *p);
   static void deleteArray__slice(void *p);
   static void destruct__slice(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::_slice*)
   {
      ::_slice *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::_slice >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("_slice", ::_slice::Class_Version(), "slice_struct.h", 25,
                  typeid(::_slice), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::_slice::Dictionary, isa_proxy, 4,
                  sizeof(::_slice) );
      instance.SetNew(&new__slice);
      instance.SetNewArray(&newArray__slice);
      instance.SetDelete(&delete__slice);
      instance.SetDeleteArray(&deleteArray__slice);
      instance.SetDestructor(&destruct__slice);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::_slice*)
   {
      return GenerateInitInstanceLocal(static_cast<::_slice*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::_slice*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr _pfp::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *_pfp::Class_Name()
{
   return "_pfp";
}

//______________________________________________________________________________
const char *_pfp::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::_pfp*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int _pfp::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::_pfp*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *_pfp::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::_pfp*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *_pfp::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::_pfp*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr _slice::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *_slice::Class_Name()
{
   return "_slice";
}

//______________________________________________________________________________
const char *_slice::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::_slice*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int _slice::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::_slice*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *_slice::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::_slice*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *_slice::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::_slice*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void _pfp::Streamer(TBuffer &R__b)
{
   // Stream an object of class _pfp.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(_pfp::Class(),this);
   } else {
      R__b.WriteClassBuffer(_pfp::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new__pfp(void *p) {
      return  p ? new(p) ::_pfp : new ::_pfp;
   }
   static void *newArray__pfp(Long_t nElements, void *p) {
      return p ? new(p) ::_pfp[nElements] : new ::_pfp[nElements];
   }
   // Wrapper around operator delete
   static void delete__pfp(void *p) {
      delete (static_cast<::_pfp*>(p));
   }
   static void deleteArray__pfp(void *p) {
      delete [] (static_cast<::_pfp*>(p));
   }
   static void destruct__pfp(void *p) {
      typedef ::_pfp current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::_pfp

//______________________________________________________________________________
void _slice::Streamer(TBuffer &R__b)
{
   // Stream an object of class _slice.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(_slice::Class(),this);
   } else {
      R__b.WriteClassBuffer(_slice::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new__slice(void *p) {
      return  p ? new(p) ::_slice : new ::_slice;
   }
   static void *newArray__slice(Long_t nElements, void *p) {
      return p ? new(p) ::_slice[nElements] : new ::_slice[nElements];
   }
   // Wrapper around operator delete
   static void delete__slice(void *p) {
      delete (static_cast<::_slice*>(p));
   }
   static void deleteArray__slice(void *p) {
      delete [] (static_cast<::_slice*>(p));
   }
   static void destruct__slice(void *p) {
      typedef ::_slice current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::_slice

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 423,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlE_slicegR_Dictionary();
   static void vectorlE_slicegR_TClassManip(TClass*);
   static void *new_vectorlE_slicegR(void *p = nullptr);
   static void *newArray_vectorlE_slicegR(Long_t size, void *p);
   static void delete_vectorlE_slicegR(void *p);
   static void deleteArray_vectorlE_slicegR(void *p);
   static void destruct_vectorlE_slicegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<_slice>*)
   {
      vector<_slice> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<_slice>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<_slice>", -2, "vector", 423,
                  typeid(vector<_slice>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlE_slicegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<_slice>) );
      instance.SetNew(&new_vectorlE_slicegR);
      instance.SetNewArray(&newArray_vectorlE_slicegR);
      instance.SetDelete(&delete_vectorlE_slicegR);
      instance.SetDeleteArray(&deleteArray_vectorlE_slicegR);
      instance.SetDestructor(&destruct_vectorlE_slicegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<_slice> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<_slice>","std::vector<_slice, std::allocator<_slice> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<_slice>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlE_slicegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<_slice>*>(nullptr))->GetClass();
      vectorlE_slicegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlE_slicegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlE_slicegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<_slice> : new vector<_slice>;
   }
   static void *newArray_vectorlE_slicegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<_slice>[nElements] : new vector<_slice>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlE_slicegR(void *p) {
      delete (static_cast<vector<_slice>*>(p));
   }
   static void deleteArray_vectorlE_slicegR(void *p) {
      delete [] (static_cast<vector<_slice>*>(p));
   }
   static void destruct_vectorlE_slicegR(void *p) {
      typedef vector<_slice> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<_slice>

namespace ROOT {
   static TClass *vectorlE_pfpgR_Dictionary();
   static void vectorlE_pfpgR_TClassManip(TClass*);
   static void *new_vectorlE_pfpgR(void *p = nullptr);
   static void *newArray_vectorlE_pfpgR(Long_t size, void *p);
   static void delete_vectorlE_pfpgR(void *p);
   static void deleteArray_vectorlE_pfpgR(void *p);
   static void destruct_vectorlE_pfpgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<_pfp>*)
   {
      vector<_pfp> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<_pfp>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<_pfp>", -2, "vector", 423,
                  typeid(vector<_pfp>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlE_pfpgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<_pfp>) );
      instance.SetNew(&new_vectorlE_pfpgR);
      instance.SetNewArray(&newArray_vectorlE_pfpgR);
      instance.SetDelete(&delete_vectorlE_pfpgR);
      instance.SetDeleteArray(&deleteArray_vectorlE_pfpgR);
      instance.SetDestructor(&destruct_vectorlE_pfpgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<_pfp> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<_pfp>","std::vector<_pfp, std::allocator<_pfp> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<_pfp>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlE_pfpgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<_pfp>*>(nullptr))->GetClass();
      vectorlE_pfpgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlE_pfpgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlE_pfpgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<_pfp> : new vector<_pfp>;
   }
   static void *newArray_vectorlE_pfpgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<_pfp>[nElements] : new vector<_pfp>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlE_pfpgR(void *p) {
      delete (static_cast<vector<_pfp>*>(p));
   }
   static void deleteArray_vectorlE_pfpgR(void *p) {
      delete [] (static_cast<vector<_pfp>*>(p));
   }
   static void destruct_vectorlE_pfpgR(void *p) {
      typedef vector<_pfp> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<_pfp>

namespace {
  void TriggerDictionaryInitialization_SliceStructDict_Impl() {
    static const char* headers[] = {
"slice_struct.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/include/",
"/exp/icarus/app/users/nsommagg/NicolaICARUS/PID/selection1muNp_PIDfinale/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "SliceStructDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$slice_struct.h")))  _pfp;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$slice_struct.h")))  _slice;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "SliceStructDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "slice_struct.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"_pfp", payloadCode, "@",
"_slice", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("SliceStructDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_SliceStructDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_SliceStructDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_SliceStructDict() {
  TriggerDictionaryInitialization_SliceStructDict_Impl();
}
