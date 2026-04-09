// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME PartonMatcher_dict
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
#include "PartonMatcher/PartonMatcher.h"
#include "PartonMatcher/PartonMatcher.h"

// Header files passed via #pragma extra_include
#include "PartonMatcher/PartonMatcher.h"

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_PartonMatcher(void *p = nullptr);
   static void *newArray_PartonMatcher(Long_t size, void *p);
   static void delete_PartonMatcher(void *p);
   static void deleteArray_PartonMatcher(void *p);
   static void destruct_PartonMatcher(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::PartonMatcher*)
   {
      ::PartonMatcher *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::PartonMatcher >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("PartonMatcher", ::PartonMatcher::Class_Version(), "PartonMatcher/PartonMatcher.h", 48,
                  typeid(::PartonMatcher), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::PartonMatcher::Dictionary, isa_proxy, 4,
                  sizeof(::PartonMatcher) );
      instance.SetNew(&new_PartonMatcher);
      instance.SetNewArray(&newArray_PartonMatcher);
      instance.SetDelete(&delete_PartonMatcher);
      instance.SetDeleteArray(&deleteArray_PartonMatcher);
      instance.SetDestructor(&destruct_PartonMatcher);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::PartonMatcher*)
   {
      return GenerateInitInstanceLocal(static_cast<::PartonMatcher*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::PartonMatcher*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr PartonMatcher::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *PartonMatcher::Class_Name()
{
   return "PartonMatcher";
}

//______________________________________________________________________________
const char *PartonMatcher::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PartonMatcher*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int PartonMatcher::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::PartonMatcher*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *PartonMatcher::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PartonMatcher*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *PartonMatcher::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::PartonMatcher*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void PartonMatcher::Streamer(TBuffer &R__b)
{
   // Stream an object of class PartonMatcher.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(PartonMatcher::Class(),this);
   } else {
      R__b.WriteClassBuffer(PartonMatcher::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_PartonMatcher(void *p) {
      return  p ? new(p) ::PartonMatcher : new ::PartonMatcher;
   }
   static void *newArray_PartonMatcher(Long_t nElements, void *p) {
      return p ? new(p) ::PartonMatcher[nElements] : new ::PartonMatcher[nElements];
   }
   // Wrapper around operator delete
   static void delete_PartonMatcher(void *p) {
      delete (static_cast<::PartonMatcher*>(p));
   }
   static void deleteArray_PartonMatcher(void *p) {
      delete [] (static_cast<::PartonMatcher*>(p));
   }
   static void destruct_PartonMatcher(void *p) {
      typedef ::PartonMatcher current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::PartonMatcher

namespace ROOT {
   // Registration Schema evolution read functions
   int RecordReadRules_libPartonMatcher() {
      return 0;
   }
   static int _R__UNIQUE_DICT_(ReadRules_libPartonMatcher) = RecordReadRules_libPartonMatcher();R__UseDummy(_R__UNIQUE_DICT_(ReadRules_libPartonMatcher));
} // namespace ROOT
namespace {
  void TriggerDictionaryInitialization_libPartonMatcher_Impl() {
    static const char* headers[] = {
"PartonMatcher/PartonMatcher.h",
nullptr
    };
    static const char* includePaths[] = {
"/data/leet/FFTutorial/atlastohyper",
"/data/leet/FFTutorial/install/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/0.6/StatAnalysis/0.6.3/InstallArea/x86_64-el9-gcc13-opt/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/0.7/StatAnalysis/0.7.2/InstallArea/x86_64-el9-gcc14-opt/include",
"/data/leet/FFTutorial/atlastohyper",
"/data/leet/FFTutorial/install/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/0.6/StatAnalysis/0.6.3/InstallArea/x86_64-el9-gcc13-opt/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/0.7/StatAnalysis/0.7.2/InstallArea/x86_64-el9-gcc14-opt/include",
"/cvmfs/atlas.cern.ch/repo/sw/software/0.6/StatAnalysis/0.6.3/InstallArea/x86_64-el9-gcc13-opt/include/",
"/data/leet/FFTutorial/atlastohyper/build_custom/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libPartonMatcher dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$PartonMatcher/PartonMatcher.h")))  PartonMatcher;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libPartonMatcher dictionary payload"

#ifndef BOOST_BIND_GLOBAL_PLACEHOLDERS
  #define BOOST_BIND_GLOBAL_PLACEHOLDERS 1
#endif
#ifndef BOOST_ALLOW_DEPRECATED_HEADERS
  #define BOOST_ALLOW_DEPRECATED_HEADERS 1
#endif
#ifndef ONNXRUNTIME_AVAILABLE
  #define ONNXRUNTIME_AVAILABLE 1
#endif
#ifndef BOOST_BIND_GLOBAL_PLACEHOLDERS
  #define BOOST_BIND_GLOBAL_PLACEHOLDERS 1
#endif
#ifndef BOOST_ALLOW_DEPRECATED_HEADERS
  #define BOOST_ALLOW_DEPRECATED_HEADERS 1
#endif
#ifndef ONNXRUNTIME_AVAILABLE
  #define ONNXRUNTIME_AVAILABLE 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "PartonMatcher/PartonMatcher.h"
#include "PartonMatcher/PartonMatcher.h"

#ifdef __CINT__

#pragma extra_include "PartonMatcher/PartonMatcher.h";
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ class PartonMatcher+;

#endif

// Extra includes
#include "PartonMatcher/PartonMatcher.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"PartonMatcher", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libPartonMatcher",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libPartonMatcher_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libPartonMatcher_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libPartonMatcher() {
  TriggerDictionaryInitialization_libPartonMatcher_Impl();
}
