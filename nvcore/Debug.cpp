// This code is in the public domain -- Ignacio Castaño <castano@gmail.com>
#include "nvcore.h"

/// Shows a message through the message handler.
void nvDebugPrint(const char *msg, ...)
{
	va_list arg;
	va_start(arg, msg);
	vprintf(msg, arg);
	va_end(arg);
}
