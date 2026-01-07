// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME daughtersInfoDict
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
#include "daughtersInfo.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_daughtersInfo(void *p = nullptr);
   static void *newArray_daughtersInfo(Long_t size, void *p);
   static void delete_daughtersInfo(void *p);
   static void deleteArray_daughtersInfo(void *p);
   static void destruct_daughtersInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::daughtersInfo*)
   {
      ::daughtersInfo *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::daughtersInfo >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("daughtersInfo", ::daughtersInfo::Class_Version(), "daughtersInfo.h", 7,
                  typeid(::daughtersInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::daughtersInfo::Dictionary, isa_proxy, 4,
                  sizeof(::daughtersInfo) );
      instance.SetNew(&new_daughtersInfo);
      instance.SetNewArray(&newArray_daughtersInfo);
      instance.SetDelete(&delete_daughtersInfo);
      instance.SetDeleteArray(&deleteArray_daughtersInfo);
      instance.SetDestructor(&destruct_daughtersInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::daughtersInfo*)
   {
      return GenerateInitInstanceLocal((::daughtersInfo*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::daughtersInfo*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr daughtersInfo::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *daughtersInfo::Class_Name()
{
   return "daughtersInfo";
}

//______________________________________________________________________________
const char *daughtersInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::daughtersInfo*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int daughtersInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::daughtersInfo*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *daughtersInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::daughtersInfo*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *daughtersInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::daughtersInfo*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void daughtersInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class daughtersInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(daughtersInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(daughtersInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_daughtersInfo(void *p) {
      return  p ? new(p) ::daughtersInfo : new ::daughtersInfo;
   }
   static void *newArray_daughtersInfo(Long_t nElements, void *p) {
      return p ? new(p) ::daughtersInfo[nElements] : new ::daughtersInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_daughtersInfo(void *p) {
      delete ((::daughtersInfo*)p);
   }
   static void deleteArray_daughtersInfo(void *p) {
      delete [] ((::daughtersInfo*)p);
   }
   static void destruct_daughtersInfo(void *p) {
      typedef ::daughtersInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::daughtersInfo

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
   static TClass *vectorlEdaughtersInfogR_Dictionary();
   static void vectorlEdaughtersInfogR_TClassManip(TClass*);
   static void *new_vectorlEdaughtersInfogR(void *p = nullptr);
   static void *newArray_vectorlEdaughtersInfogR(Long_t size, void *p);
   static void delete_vectorlEdaughtersInfogR(void *p);
   static void deleteArray_vectorlEdaughtersInfogR(void *p);
   static void destruct_vectorlEdaughtersInfogR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<daughtersInfo>*)
   {
      vector<daughtersInfo> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<daughtersInfo>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<daughtersInfo>", -2, "vector", 386,
                  typeid(vector<daughtersInfo>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdaughtersInfogR_Dictionary, isa_proxy, 4,
                  sizeof(vector<daughtersInfo>) );
      instance.SetNew(&new_vectorlEdaughtersInfogR);
      instance.SetNewArray(&newArray_vectorlEdaughtersInfogR);
      instance.SetDelete(&delete_vectorlEdaughtersInfogR);
      instance.SetDeleteArray(&deleteArray_vectorlEdaughtersInfogR);
      instance.SetDestructor(&destruct_vectorlEdaughtersInfogR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<daughtersInfo> >()));

      ::ROOT::AddClassAlternate("vector<daughtersInfo>","std::vector<daughtersInfo, std::allocator<daughtersInfo> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<daughtersInfo>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdaughtersInfogR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<daughtersInfo>*)nullptr)->GetClass();
      vectorlEdaughtersInfogR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdaughtersInfogR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdaughtersInfogR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<daughtersInfo> : new vector<daughtersInfo>;
   }
   static void *newArray_vectorlEdaughtersInfogR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<daughtersInfo>[nElements] : new vector<daughtersInfo>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdaughtersInfogR(void *p) {
      delete ((vector<daughtersInfo>*)p);
   }
   static void deleteArray_vectorlEdaughtersInfogR(void *p) {
      delete [] ((vector<daughtersInfo>*)p);
   }
   static void destruct_vectorlEdaughtersInfogR(void *p) {
      typedef vector<daughtersInfo> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<daughtersInfo>

namespace {
  void TriggerDictionaryInitialization_daughtersInfoDict_Impl() {
    static const char* headers[] = {
"daughtersInfo.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_26_06b/Linux64bit+3.10-2.17-e20-p3913-prof/include/",
"/storage/gpfs_data/icarus/local/users/sommaggio/simul_z/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "daughtersInfoDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate(R"ATTRDUMP(ROOT macro per la dictionary)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$daughtersInfo.h")))  daughtersInfo;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "daughtersInfoDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "daughtersInfo.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"daughtersInfo", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("daughtersInfoDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_daughtersInfoDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_daughtersInfoDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_daughtersInfoDict() {
  TriggerDictionaryInitialization_daughtersInfoDict_Impl();
}
