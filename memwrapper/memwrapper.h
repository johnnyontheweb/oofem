// memwrapper.h
#include <memory>
#include <string>

#pragma once
#pragma unmanaged

#if defined(LIB_EXPORT)
#define DECLSPEC_CLASS __declspec(dllexport)
#else
#define DECLSPEC_CLASS __declspec(dllimport)
#endif

//using namespace System;
#ifndef MEMSTR
namespace oofem {
#endif MEMSTR
	class MemStrPrivate;

	class __declspec(dllexport) MemWrapper
	{

		static std::unique_ptr<MemStrPrivate> s_instance;
		static std::string oofemPID;

	public:
		// MemWrapper();
		static MemStrPrivate* instance();// const std::auto_ptr<char*> &pid);

		static FILE* getInstance(const char * key, const std::string &pid);
		//~MemWrapper();
	};
#ifndef MEMSTR
}

#endif MEMSTR