// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME DataStructDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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
#include "data_struct.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_reco_track(void *p = nullptr);
   static void *newArray_reco_track(Long_t size, void *p);
   static void delete_reco_track(void *p);
   static void deleteArray_reco_track(void *p);
   static void destruct_reco_track(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::reco_track*)
   {
      ::reco_track *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::reco_track >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("reco_track", ::reco_track::Class_Version(), "data_struct.h", 10,
                  typeid(::reco_track), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::reco_track::Dictionary, isa_proxy, 4,
                  sizeof(::reco_track) );
      instance.SetNew(&new_reco_track);
      instance.SetNewArray(&newArray_reco_track);
      instance.SetDelete(&delete_reco_track);
      instance.SetDeleteArray(&deleteArray_reco_track);
      instance.SetDestructor(&destruct_reco_track);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::reco_track*)
   {
      return GenerateInitInstanceLocal((::reco_track*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::reco_track*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_true_track(void *p = nullptr);
   static void *newArray_true_track(Long_t size, void *p);
   static void delete_true_track(void *p);
   static void deleteArray_true_track(void *p);
   static void destruct_true_track(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::true_track*)
   {
      ::true_track *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::true_track >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("true_track", ::true_track::Class_Version(), "data_struct.h", 48,
                  typeid(::true_track), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::true_track::Dictionary, isa_proxy, 4,
                  sizeof(::true_track) );
      instance.SetNew(&new_true_track);
      instance.SetNewArray(&newArray_true_track);
      instance.SetDelete(&delete_true_track);
      instance.SetDeleteArray(&deleteArray_true_track);
      instance.SetDestructor(&destruct_true_track);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::true_track*)
   {
      return GenerateInitInstanceLocal((::true_track*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::true_track*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_track(void *p = nullptr);
   static void *newArray_track(Long_t size, void *p);
   static void delete_track(void *p);
   static void deleteArray_track(void *p);
   static void destruct_track(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::track*)
   {
      ::track *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::track >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("track", ::track::Class_Version(), "data_struct.h", 63,
                  typeid(::track), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::track::Dictionary, isa_proxy, 4,
                  sizeof(::track) );
      instance.SetNew(&new_track);
      instance.SetNewArray(&newArray_track);
      instance.SetDelete(&delete_track);
      instance.SetDeleteArray(&deleteArray_track);
      instance.SetDestructor(&destruct_track);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::track*)
   {
      return GenerateInitInstanceLocal((::track*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::track*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_RecoSlice(void *p = nullptr);
   static void *newArray_RecoSlice(Long_t size, void *p);
   static void delete_RecoSlice(void *p);
   static void deleteArray_RecoSlice(void *p);
   static void destruct_RecoSlice(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RecoSlice*)
   {
      ::RecoSlice *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RecoSlice >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("RecoSlice", ::RecoSlice::Class_Version(), "data_struct.h", 77,
                  typeid(::RecoSlice), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RecoSlice::Dictionary, isa_proxy, 4,
                  sizeof(::RecoSlice) );
      instance.SetNew(&new_RecoSlice);
      instance.SetNewArray(&newArray_RecoSlice);
      instance.SetDelete(&delete_RecoSlice);
      instance.SetDeleteArray(&deleteArray_RecoSlice);
      instance.SetDestructor(&destruct_RecoSlice);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RecoSlice*)
   {
      return GenerateInitInstanceLocal((::RecoSlice*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RecoSlice*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr reco_track::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *reco_track::Class_Name()
{
   return "reco_track";
}

//______________________________________________________________________________
const char *reco_track::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::reco_track*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int reco_track::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::reco_track*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *reco_track::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::reco_track*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *reco_track::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::reco_track*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr true_track::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *true_track::Class_Name()
{
   return "true_track";
}

//______________________________________________________________________________
const char *true_track::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::true_track*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int true_track::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::true_track*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *true_track::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::true_track*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *true_track::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::true_track*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr track::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *track::Class_Name()
{
   return "track";
}

//______________________________________________________________________________
const char *track::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::track*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int track::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::track*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *track::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::track*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *track::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::track*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr RecoSlice::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *RecoSlice::Class_Name()
{
   return "RecoSlice";
}

//______________________________________________________________________________
const char *RecoSlice::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RecoSlice*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int RecoSlice::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RecoSlice*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RecoSlice::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RecoSlice*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RecoSlice::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RecoSlice*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void reco_track::Streamer(TBuffer &R__b)
{
   // Stream an object of class reco_track.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(reco_track::Class(),this);
   } else {
      R__b.WriteClassBuffer(reco_track::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_reco_track(void *p) {
      return  p ? new(p) ::reco_track : new ::reco_track;
   }
   static void *newArray_reco_track(Long_t nElements, void *p) {
      return p ? new(p) ::reco_track[nElements] : new ::reco_track[nElements];
   }
   // Wrapper around operator delete
   static void delete_reco_track(void *p) {
      delete ((::reco_track*)p);
   }
   static void deleteArray_reco_track(void *p) {
      delete [] ((::reco_track*)p);
   }
   static void destruct_reco_track(void *p) {
      typedef ::reco_track current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::reco_track

//______________________________________________________________________________
void true_track::Streamer(TBuffer &R__b)
{
   // Stream an object of class true_track.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(true_track::Class(),this);
   } else {
      R__b.WriteClassBuffer(true_track::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_true_track(void *p) {
      return  p ? new(p) ::true_track : new ::true_track;
   }
   static void *newArray_true_track(Long_t nElements, void *p) {
      return p ? new(p) ::true_track[nElements] : new ::true_track[nElements];
   }
   // Wrapper around operator delete
   static void delete_true_track(void *p) {
      delete ((::true_track*)p);
   }
   static void deleteArray_true_track(void *p) {
      delete [] ((::true_track*)p);
   }
   static void destruct_true_track(void *p) {
      typedef ::true_track current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::true_track

//______________________________________________________________________________
void track::Streamer(TBuffer &R__b)
{
   // Stream an object of class track.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(track::Class(),this);
   } else {
      R__b.WriteClassBuffer(track::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_track(void *p) {
      return  p ? new(p) ::track : new ::track;
   }
   static void *newArray_track(Long_t nElements, void *p) {
      return p ? new(p) ::track[nElements] : new ::track[nElements];
   }
   // Wrapper around operator delete
   static void delete_track(void *p) {
      delete ((::track*)p);
   }
   static void deleteArray_track(void *p) {
      delete [] ((::track*)p);
   }
   static void destruct_track(void *p) {
      typedef ::track current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::track

//______________________________________________________________________________
void RecoSlice::Streamer(TBuffer &R__b)
{
   // Stream an object of class RecoSlice.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(RecoSlice::Class(),this);
   } else {
      R__b.WriteClassBuffer(RecoSlice::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RecoSlice(void *p) {
      return  p ? new(p) ::RecoSlice : new ::RecoSlice;
   }
   static void *newArray_RecoSlice(Long_t nElements, void *p) {
      return p ? new(p) ::RecoSlice[nElements] : new ::RecoSlice[nElements];
   }
   // Wrapper around operator delete
   static void delete_RecoSlice(void *p) {
      delete ((::RecoSlice*)p);
   }
   static void deleteArray_RecoSlice(void *p) {
      delete [] ((::RecoSlice*)p);
   }
   static void destruct_RecoSlice(void *p) {
      typedef ::RecoSlice current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::RecoSlice

namespace ROOT {
   static TClass *vectorlEtrue_trackgR_Dictionary();
   static void vectorlEtrue_trackgR_TClassManip(TClass*);
   static void *new_vectorlEtrue_trackgR(void *p = nullptr);
   static void *newArray_vectorlEtrue_trackgR(Long_t size, void *p);
   static void delete_vectorlEtrue_trackgR(void *p);
   static void deleteArray_vectorlEtrue_trackgR(void *p);
   static void destruct_vectorlEtrue_trackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<true_track>*)
   {
      vector<true_track> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<true_track>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<true_track>", -2, "vector", 386,
                  typeid(vector<true_track>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrue_trackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<true_track>) );
      instance.SetNew(&new_vectorlEtrue_trackgR);
      instance.SetNewArray(&newArray_vectorlEtrue_trackgR);
      instance.SetDelete(&delete_vectorlEtrue_trackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrue_trackgR);
      instance.SetDestructor(&destruct_vectorlEtrue_trackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<true_track> >()));

      ::ROOT::AddClassAlternate("vector<true_track>","std::vector<true_track, std::allocator<true_track> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<true_track>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrue_trackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<true_track>*)nullptr)->GetClass();
      vectorlEtrue_trackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrue_trackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrue_trackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<true_track> : new vector<true_track>;
   }
   static void *newArray_vectorlEtrue_trackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<true_track>[nElements] : new vector<true_track>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrue_trackgR(void *p) {
      delete ((vector<true_track>*)p);
   }
   static void deleteArray_vectorlEtrue_trackgR(void *p) {
      delete [] ((vector<true_track>*)p);
   }
   static void destruct_vectorlEtrue_trackgR(void *p) {
      typedef vector<true_track> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<true_track>

namespace ROOT {
   static TClass *vectorlEtrackgR_Dictionary();
   static void vectorlEtrackgR_TClassManip(TClass*);
   static void *new_vectorlEtrackgR(void *p = nullptr);
   static void *newArray_vectorlEtrackgR(Long_t size, void *p);
   static void delete_vectorlEtrackgR(void *p);
   static void deleteArray_vectorlEtrackgR(void *p);
   static void destruct_vectorlEtrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<track>*)
   {
      vector<track> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<track>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<track>", -2, "vector", 386,
                  typeid(vector<track>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEtrackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<track>) );
      instance.SetNew(&new_vectorlEtrackgR);
      instance.SetNewArray(&newArray_vectorlEtrackgR);
      instance.SetDelete(&delete_vectorlEtrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEtrackgR);
      instance.SetDestructor(&destruct_vectorlEtrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<track> >()));

      ::ROOT::AddClassAlternate("vector<track>","std::vector<track, std::allocator<track> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<track>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEtrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<track>*)nullptr)->GetClass();
      vectorlEtrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEtrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEtrackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<track> : new vector<track>;
   }
   static void *newArray_vectorlEtrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<track>[nElements] : new vector<track>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEtrackgR(void *p) {
      delete ((vector<track>*)p);
   }
   static void deleteArray_vectorlEtrackgR(void *p) {
      delete [] ((vector<track>*)p);
   }
   static void destruct_vectorlEtrackgR(void *p) {
      typedef vector<track> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<track>

namespace ROOT {
   static TClass *vectorlEreco_trackgR_Dictionary();
   static void vectorlEreco_trackgR_TClassManip(TClass*);
   static void *new_vectorlEreco_trackgR(void *p = nullptr);
   static void *newArray_vectorlEreco_trackgR(Long_t size, void *p);
   static void delete_vectorlEreco_trackgR(void *p);
   static void deleteArray_vectorlEreco_trackgR(void *p);
   static void destruct_vectorlEreco_trackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<reco_track>*)
   {
      vector<reco_track> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<reco_track>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<reco_track>", -2, "vector", 386,
                  typeid(vector<reco_track>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEreco_trackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<reco_track>) );
      instance.SetNew(&new_vectorlEreco_trackgR);
      instance.SetNewArray(&newArray_vectorlEreco_trackgR);
      instance.SetDelete(&delete_vectorlEreco_trackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEreco_trackgR);
      instance.SetDestructor(&destruct_vectorlEreco_trackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<reco_track> >()));

      ::ROOT::AddClassAlternate("vector<reco_track>","std::vector<reco_track, std::allocator<reco_track> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<reco_track>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEreco_trackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<reco_track>*)nullptr)->GetClass();
      vectorlEreco_trackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEreco_trackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEreco_trackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<reco_track> : new vector<reco_track>;
   }
   static void *newArray_vectorlEreco_trackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<reco_track>[nElements] : new vector<reco_track>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEreco_trackgR(void *p) {
      delete ((vector<reco_track>*)p);
   }
   static void deleteArray_vectorlEreco_trackgR(void *p) {
      delete [] ((vector<reco_track>*)p);
   }
   static void destruct_vectorlEreco_trackgR(void *p) {
      typedef vector<reco_track> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<reco_track>

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
         instance("vector<double>", -2, "vector", 386,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)nullptr)->GetClass();
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
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlERecoSlicegR_Dictionary();
   static void vectorlERecoSlicegR_TClassManip(TClass*);
   static void *new_vectorlERecoSlicegR(void *p = nullptr);
   static void *newArray_vectorlERecoSlicegR(Long_t size, void *p);
   static void delete_vectorlERecoSlicegR(void *p);
   static void deleteArray_vectorlERecoSlicegR(void *p);
   static void destruct_vectorlERecoSlicegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RecoSlice>*)
   {
      vector<RecoSlice> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RecoSlice>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RecoSlice>", -2, "vector", 386,
                  typeid(vector<RecoSlice>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERecoSlicegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<RecoSlice>) );
      instance.SetNew(&new_vectorlERecoSlicegR);
      instance.SetNewArray(&newArray_vectorlERecoSlicegR);
      instance.SetDelete(&delete_vectorlERecoSlicegR);
      instance.SetDeleteArray(&deleteArray_vectorlERecoSlicegR);
      instance.SetDestructor(&destruct_vectorlERecoSlicegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RecoSlice> >()));

      ::ROOT::AddClassAlternate("vector<RecoSlice>","std::vector<RecoSlice, std::allocator<RecoSlice> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RecoSlice>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERecoSlicegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RecoSlice>*)nullptr)->GetClass();
      vectorlERecoSlicegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERecoSlicegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERecoSlicegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RecoSlice> : new vector<RecoSlice>;
   }
   static void *newArray_vectorlERecoSlicegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RecoSlice>[nElements] : new vector<RecoSlice>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERecoSlicegR(void *p) {
      delete ((vector<RecoSlice>*)p);
   }
   static void deleteArray_vectorlERecoSlicegR(void *p) {
      delete [] ((vector<RecoSlice>*)p);
   }
   static void destruct_vectorlERecoSlicegR(void *p) {
      typedef vector<RecoSlice> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RecoSlice>

namespace {
  void TriggerDictionaryInitialization_DataStructDict_Impl() {
    static const char* headers[] = {
"data_struct.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/include/",
"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/PID/selection_new_PID/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "DataStructDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$data_struct.h")))  track;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$data_struct.h")))  RecoSlice;
class __attribute__((annotate("$clingAutoload$data_struct.h")))  true_track;
class __attribute__((annotate("$clingAutoload$data_struct.h")))  reco_track;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "DataStructDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "data_struct.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"RecoSlice", payloadCode, "@",
"reco_track", payloadCode, "@",
"track", payloadCode, "@",
"true_track", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DataStructDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DataStructDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DataStructDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DataStructDict() {
  TriggerDictionaryInitialization_DataStructDict_Impl();
}
