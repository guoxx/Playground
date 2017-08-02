#include <Windows.h>
#include "RenderDoc.h"
#include "renderdoc_app.h"

static RENDERDOC_API_1_1_1* s_renderdoc{nullptr};

void RenderDoc::initialize(bool debug)
{
	s_renderdoc = new RENDERDOC_API_1_1_1{};

	HMODULE renderdoc = LoadLibrary(TEXT("renderdoc.dll"));
	pRENDERDOC_GetAPI renderdocGetApiFunc = (pRENDERDOC_GetAPI)GetProcAddress(renderdoc, "RENDERDOC_GetAPI");
	renderdocGetApiFunc(eRENDERDOC_API_Version_1_1_1, reinterpret_cast<void**>(&s_renderdoc));

	s_renderdoc->SetLogFilePathTemplate("renderdoc/playground");
	if (debug)
	{
		s_renderdoc->SetCaptureOptionU32(eRENDERDOC_Option_APIValidation, 1);
		s_renderdoc->SetCaptureOptionU32(eRENDERDOC_Option_DebugOutputMute, 0);
	}

    s_renderdoc->SetCaptureOptionU32(eRENDERDOC_Option_CaptureCallstacks, 1);

	RENDERDOC_InputButton k = eRENDERDOC_Key_PrtScrn;
	s_renderdoc->SetCaptureKeys(&k, 1);
}

void RenderDoc::finalize()
{
	delete s_renderdoc;
	s_renderdoc = nullptr;
}

