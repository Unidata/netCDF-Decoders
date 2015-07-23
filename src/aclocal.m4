dnl Set ulog parameters
dnl
AC_DEFUN([UL_ULOG], [dnl
    AC_MSG_CHECKING(ulog defines)
    CPP_ULOG=
    case `uname -s` in
    OSF1|sn1036|Linux)
	CPP_ULOG=-DNO_REPLACE_SYSLOG
	;;
    esac
    if test -w /dev/log; then
	if file /dev/log | grep socket >/dev/null; then
	    CPP_ULOG="${CPP_ULOG} -DLOGNAME_ISSOCK"
    	fi
    elif test -w /dev/conslog; then
        CPP_ULOG="${CPP_ULOG} -D_DEV_CONSLOG"
    fi
    AC_SUBST(CPP_ULOG)
    AC_MSG_RESULT($CPP_ULOG)
])
dnl Set Berkeley socket references.
dnl
