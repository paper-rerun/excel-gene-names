# Drop the current process to low priority, so the run does not get in the way of the desktop
import os
import sys

BELOW_NORMAL_PRIORITY_CLASS = 0x00004000


def SetLowPriority():
	if sys.platform == "win32":
		import ctypes
		# -1 is the pseudo handle of the current process, child processes inherit the class
		ctypes.windll.kernel32.SetPriorityClass(ctypes.c_void_p(-1), BELOW_NORMAL_PRIORITY_CLASS)
	else:
		os.nice(19)


def GetPriorityClass():
	if sys.platform == "win32":
		import ctypes
		return ctypes.windll.kernel32.GetPriorityClass(ctypes.c_void_p(-1))
	return os.nice(0)
