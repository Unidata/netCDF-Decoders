/*
 *   Copyright 1993, University Corporation for Atmospheric Research
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: ulog.c,v 1.3 2003/02/27 21:47:47 rkambic Exp $ */

/* 
 * Utility Functions to implement consistant logging mechanisms.
 * These provide interfaces to the syslogd(8) system if available.
 * 'syslog' code based on "@(#)syslog.c	5.20 (Berkeley) 1/19/89"
 * Copyright (c) 1983, 1988 Regents of the University of California.
 */

/*
 * If a file /dev/conslog exists, configure defines _DEV_CONSLOG.
 * We use this macro to handle SVR4 streams based logging.
 */
#ifndef _DEV_CONSLOG
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE
#endif  /* !_POSIX_SOURCE */
#endif /* !_DEV_CONSLOG */

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#if defined(__STDC__) || defined(_AIX)
#include <stdarg.h>
#define STDC_ARGS 1
#else
#include <varargs.h>
#endif
#include <errno.h>
#include <time.h>
#include <signal.h>
#ifndef NO_WAITPID
#include <sys/wait.h>
#endif
#include "ulog.h"

#if defined(__CENTERLINE__) && defined(sun) && defined(__STDC__)
/* Workaround for ObjectCenter 1.1 stdargs problem, Sun OS 4.1.x */
#undef __builtin_va_arg_incr
#define __builtin_va_arg_incr(list) (((list) += 1) -1)
#undef va_start
#define va_start(list, arg) list = (char *)&arg + sizeof(arg)
#endif

#if defined(NO_STRERROR) && !defined(lint) /* recividist unie */
char *
strerror(errnum)
int errnum;
{
    extern int sys_nerr;
    extern char *sys_errlist[];

    if(errnum < 0 || errnum >= sys_nerr) return NULL;
    /* else */
    return sys_errlist[errnum];
}
#elif defined(sun) /* Sun OS 4.x */
extern char *strerror(/* int */);
#endif /* NO_STRERROR */


#if (LOG_NFACILITIES > 0)
/*
 * BSD 4.3 style, the logger attempts to open a connection to
 * a named fifo or an AF_UNIX socket called ULOGNAME. 
 */
#ifdef _DEV_CONSLOG
#define	ULOGNAME	"/dev/conslog"
#else
#define	ULOGNAME	"/dev/log"
#endif /* !_DEV_CONSLOG */
#endif

/*
 * BSD 4.2 style, the logger trys the AF_INET well known service
 * "syslog" at host LOGHOST. 
 */
#define	LOGHOST	"127.0.0.1"


/*
 * If all else fails, it attempts to write directly to CONSOLE.
 */ 
#define	CONSOLE	"/dev/console"


static const char	*logFilename = NULL;	/* Logfile name, == NULL => use syslogd */
static int	logFd = -1;		/* fd for log */
static int	logOptions = 0;		/* option bits, set by openulog() */
static char	logIdent[LOG_IDENT_LEN +1] = "ulog";	/* string to tag the entry with */
static int	logFacility = LOG_USER;		/* default facility code */
static int	logMask =  LOG_UPTO(LOG_DEBUG);		/* default mask */


/*
 * Close the log connection
 */
int
closeulog(void)
{
	int ret = 0;

	if(logFilename != NULL && *logFilename == '-') /* special case for stderr */
		logFilename = NULL;
	else if (logFd >= 0)
	{
		ret = close(logFd);
	}
	logFd = -1;

	return ret;
}


/*
 * Set the ident field.
 * Use this instead of closeulog() followed by openulog()
 * to change the ident in a child proc.
 */
void
setulogident(const char *ident)
{
	if (ident != NULL)
		strncpy(logIdent, ident, LOG_IDENT_LEN);
}


/* 
 * Like Berkeley 'openlog(3)' - initialize syslog system.
 * But: 
 *  Always attempts to open the log connection if there is none,
 *  even if the LOG_NDELAY option bit is cleared.
 *  Returns the descriptor (socket) of the logger or -1 on system call error.
 * N.B. multiple calls without an intervening 'closeulog()' simply reinitialize
 * ident, options, and facility.
 * The data referred to by 'ident' and 'filename' should have process lifetime.
 */
int
openulog(
	const char *ident,
	int options,
	int facility,
	const char *logfilename)
{
	setulogident(ident);

	logOptions = options;

	if (facility != 0 && (facility &~ LOG_FACMASK) == 0)
		logFacility = facility;

	if (logFd != -1)
		(void) closeulog();
	
	if(logfilename != NULL && *logfilename != 0)
	{
			/* "-" means use stderr */
		if(logfilename[0] == '-' && logfilename[1] == 0)
		{
			logFd = fileno(stderr);
			logFilename = "-";
		}
		else
		{
			logFd = open(logfilename, (O_WRONLY|O_CREAT|O_APPEND), 0664);
			if(logFd != -1)
				logFilename = logfilename; /* N.B. Not a copy */
		}
		/* if the open fail, fall thru and use syslogd */
	}

	if (logFd == -1)
	{
		/*
		 * On a given OS, the ULOGNAME is always
		 * a fifo, a STREAMS device,  or a unix domain socket.
		 * So, figure it out ahead of time.
		 */
#if defined(ULOGNAME)
#if LOGNAME_ISSOCK
		logFd = usopen(ULOGNAME);	/* unix domain socket open */
#else
		logFd = open(ULOGNAME, O_WRONLY);	/* System V */
#endif /*LOGNAME_ISSOCK*/
#else
		logFd = udpopen(LOGHOST, "syslog"); /* UDP socket open */
#endif /*ULOGNAME*/
	}

	if (logFd == -1)
		return -1;
	/* else */

#ifdef FD_CLOEXEC /* only try to do this on systems that support it */
	/* set descriptor to "close on exec" */
	if(logFilename == NULL || logFd != 2) /* special case for stderr */
	{
	    if (fcntl(logFd, F_SETFD, FD_CLOEXEC) == -1)
	    {
		    (void) close(logFd);
		    return -1;
	    }
	    /* else */
	}
#endif

	return logFd;
}


#ifndef LOG_PRI
#define LOG_PRI(p)	((p) & LOG_PRIMASK)	/* extract priority */
#endif
#ifndef LOG_FAC
#define LOG_FAC(p)	(((p) & LOG_FACMASK) >> 3)	/* facility of pri */
#endif


#ifdef _DEV_CONSLOG
#include <stropts.h> /* putmsg, struct strbuf */
#include <sys/strlog.h> /* struct log_ctl */
#endif 
/*
 * analogous to vsyslog()
 */
#ifdef STDC_ARGS
static int
vulog(int pri, const char *fmt, va_list args)
{
#else
int
vulog(pri, fmt, args)
	int pri;
	char *fmt;
    va_list args;
{
#endif /* STDC_ARGS */
	char tbuf[2048];
	char *cp;
	time_t now;
	size_t cnt;
	pid_t pid;

#if LOG_NFACILITIES > 0
	/* sanity check */
	if ((unsigned int)LOG_FAC(pri) >= LOG_NFACILITIES)
		return 0;
#endif

#if MASKDEBUG
	printf("                   pri 0x%08x\n", pri);
	printf("          LOG_PRI(pri) 0x%08x\n", LOG_PRI(pri));
	printf("LOG_MASK(LOG_PRI(pri)) 0x%08x\n", LOG_MASK(LOG_PRI(pri)));
	printf("               logMask 0x%08x\n", logMask);
	printf("                result 0x%08x\n",
		(LOG_MASK(LOG_PRI(pri)) & logMask));
#endif

	/* see if we should just throw out this message */
	if ((LOG_MASK(LOG_PRI(pri)) & logMask) == 0 )
		return 0;

#if LOG_NFACILITIES > 0
	/* set default facility if none specified */
	if ((pri & LOG_FACMASK) == 0)
		pri |= logFacility;
#endif

	/*
	 * If this open fails, we will eventually write the message to the console.
	 */
	if (logFd == -1)
		(void) openulog(logIdent, logOptions, logFacility, logFilename);

	/*
	 * build the message
	 */
	cp = tbuf;
#ifndef _DEV_CONSLOG
	if(logFilename == NULL) /* using syslogd */
	{
		(void) sprintf(tbuf, "<%d>", pri);
		for (; *cp != 0; ++cp)
			/*EMPTY*/;
	}
#endif /* !_DEV_CONSLOG */
		
	if(!(logOptions & LOG_NOTIME))
	{
		struct tm tm_now;
		struct tm *tmp = NULL;
		(void) time(&now);

		/* N.B.: default for this package is to use gmt */
		if(!(logOptions & LOG_LOCALTIME))
			tmp = gmtime(&now);
		/* else, also covers failure of gmtime call */
		if(!tmp)
			tmp = localtime(&now);
		tm_now = *tmp;
		/* (void) sprintf(cp, "%.15s ", ctime(&now) + 4); */
		(void) strftime(cp, 17, "%b %d %H:%M:%S ", &tm_now);

		for (; *cp != 0; ++cp)
			/*EMPTY*/;
		errno = 0; /* ultrix 4.0 mips trashes errno in ctime */
	}

	if (logIdent)
	{
		(void) strcpy(cp, logIdent);
		for (; *cp != 0; ++cp)
			/*EMPTY*/;
	}

	if (logOptions & LOG_PID)
	{
		(void) sprintf(cp, "[%d]", getpid());
		for (; *cp != 0; ++cp)
			/*EMPTY*/;
	}

	if (logIdent)
	{
		*cp++ = ':';
		*cp++ = ' ';
	}

	/* we don't do %m substitution */

	(void)vsprintf(cp, fmt, args);

	cnt = strlen(tbuf); /* cnt for write() below */
#ifdef _DEV_CONSLOG
	if(logFilename == NULL) /* using syslogd */
	{
		/* putmsg fails ERANGE when the message is too long */
		if(cnt > 126)
		{
			tbuf[127] = '\0';
			cnt = 126;
		}
	}
#endif
	if(tbuf[cnt -1] != '\n')
	{
		tbuf[cnt++] = '\n';
		tbuf[cnt] = 0;
	}
#if TBUFDIAG
	fputs(tbuf, stdout);
#endif

	/* output the message to the logger */
#ifdef _DEV_CONSLOG
	if(logFilename == NULL) /* using syslogd */
	{
		struct strbuf ctl, dat;
		struct log_ctl lc;

		lc.level = 0;
		lc.flags = (SL_NOTIFY|SL_ERROR|SL_CONSOLE|SL_WARN|SL_NOTE);
		lc.pri = pri;

		ctl.len = ctl.maxlen = sizeof(lc);
		ctl.buf = (char *)&lc;

		dat.len = dat.maxlen = cnt;
		dat.buf = tbuf;

again:
		errno = 0;
		if(putmsg(logFd, &ctl, &dat, 0) < 0)
		{
			if(errno == EAGAIN || errno == EINTR)
				goto again; /* Arrgh System V */
			(void) closeulog();
		}
	}
	else /* write to a file */
#endif /* !_DEV_CONSLOG */
	if (logFd != -1 )
	{
		if (write(logFd, tbuf, cnt) < 0)
			(void) closeulog();
	}

	if(logFd != -1 || !(logOptions & LOG_CONS) )
			return (int)cnt;
	/* else */

	/* output the message to the console */
	pid = fork();
	if (pid == -1) /* fork error */
		return -1;
	/* else */

	if (pid == 0)
	{
		/* child */
		int fd;

		(void)signal(SIGALRM, SIG_DFL);
#ifdef SIG_SETMASK /* only try to do this on systems that support it */
		{
			/* sigsetmask((long)~sigmask(SIGALRM)); */
			sigset_t set;
			(void) sigfillset(&set);
			(void) sigdelset(&set, SIGALRM);
			(void) sigprocmask(SIG_SETMASK, &set, (sigset_t *)0);
		}
#endif
		(void)alarm((unsigned int)5);
		if ((fd = open(CONSOLE, O_WRONLY, 0)) < 0)
		{
			logOptions &= (~LOG_CONS); /* don't keep trying */
			_exit(1);
		}

		/* (void)strcat(tbuf, "\r"); cnt++; */
		if(tbuf[0] == '<')
		{
			cp = strchr(tbuf, '>');
			if(cp == NULL)
				cp = tbuf;
			else
				cp++;
		}
		else
			cp = tbuf;


		(void)write(fd, cp, cnt - (size_t)(cp - tbuf));
		(void)close(fd);
		_exit(0);
	}
	/* else, parent */

	if (!(logOptions & LOG_NOWAIT))
#ifdef NO_WAITPID
		while ((cnt = wait((int *)0)) > 0 && cnt != pid) /*EMPTY*/;
#else
		(void) waitpid(pid, (int *)0, 0);	
#endif
	return 0;
}


/*
 * analogous to syslog()
 */
#ifdef STDC_ARGS
void
ulog(int pri, const char *fmt, ...)
{
#else
/*VARARGS2*/
void
ulog(pri, fmt, va_alist)
	int pri;
	char *fmt;
	va_dcl
{
#endif /* STDC_ARGS */
    va_list args;
#ifdef STDC_ARGS
	va_start(args, fmt);
#else
	va_start(args);
#endif /* !STDC_ARGS */
	(void)vulog(pri, fmt, args);
	va_end(args);
}


#ifndef NO_REPLACE_SYSLOG
/*
 * replace syslog (you may not want to do this...)
 */
#ifdef STDC_ARGS
void
syslog(int pri, const char *fmt, ...)
{
#else
/*VARARGS2*/
void
syslog(pri, fmt, va_alist)
	int pri;
	char *fmt;
	va_dcl
{
#endif /* STDC_ARGS */
    va_list args;
#ifdef STDC_ARGS
	va_start(args, fmt);
#else
	va_start(args);
#endif /* !STDC_ARGS */
	(void)vulog(pri, fmt, args);
	va_end(args);
}
#endif /* NO_REPLACE_SYSLOG */


/*
 * Set the ulog mask, return the old mask.
 * Analogous to setlogmask(), execept it actually does something.
 */
int
setulogmask(int pmask)
{
	int omask;
	
	omask = logMask;
	if (pmask != 0)
		logMask = pmask;
	return (omask);
}


/*
 * If the bit in the logMask corresponding to * pri is set,
 *    unset it.
 * Otherwise, set it.
 * Used for toggling the verbosity, eg
 * toggleulogpri(LOG_INFO);
 */
int
toggleulogpri(int pri)
{
	pri = LOG_MASK(pri);
	if(pri & logMask)
	{
		logMask &= ~pri;
		return 0;
	}
	/* else */
	logMask |= pri;
	return 1;
}


/*
 * Cycle through logging priorities:
 *  Silent, Verbose (LOG_INFO), Verbose and debug.
 */
void
rollulogpri(void)
{
	switch(logMask & (LOG_MASK(LOG_INFO)|LOG_MASK(LOG_DEBUG))) {

	case (LOG_MASK(LOG_INFO)|LOG_MASK(LOG_DEBUG)):
		unotice("Going silent");
		logMask &= ~(LOG_MASK(LOG_INFO)|LOG_MASK(LOG_DEBUG));
		break;
	case LOG_MASK(LOG_INFO):
		unotice("Adding debug");
		logMask |= LOG_MASK(LOG_DEBUG);
		break;
	case LOG_MASK(LOG_DEBUG):
		unotice("Adding verbose");
		logMask |= LOG_MASK(LOG_INFO);
		break;
	case 0:
		unotice("Going verbose");
		logMask |= LOG_MASK(LOG_INFO);
		break;
	}
}

/*
 * Get the ulog mask.
 * Useful for determining verbosity.
 */
int
getulogmask(void)
{
	return (logMask);
}


int
ulogIsVerbose(void)
{
	return (logMask & LOG_MASK(LOG_INFO));
}


int
ulogIsDebug(void)
{
	return (logMask & LOG_MASK(LOG_DEBUG));
}


/*
 * Log system call errors
 * Use where you would want to call perror(3).
 * Calling sequence is
 *	serror(format, arg1, arg2,...)
 * with zero or more args of types compatible with the associated format
 * specifiers.  For example:
 *         serror("shutting down");
 *	   serror("can't open %s", file_name);
 *         serror("process %d in state %s",pid,state);
 */
#ifdef STDC_ARGS
void
serror(const char *fmt, ...)
#else
/*VARARGS1*/
void
serror(fmt, va_alist)
     char *fmt;
     va_dcl
#endif /* STDC_ARGS */
{
    va_list args;
#ifdef STDC_ARGS
	va_start(args ,fmt);
#else
	va_start(args);
#endif /* STDC_ARGS */
	if(errno != 0)
	{
		char buf[1024];
		int errnum = errno;		/* save real errno in case we wipe it out */
		char *cp;
		(void) vsprintf(buf, fmt, args);
		for(cp = buf; *cp != 0; cp++) /*EMPTY*/;
		strcat(cp, ": %s");
		ulog(LOG_ERR, buf, strerror(errnum));
		errno = 0;
	}
	else
	{
		(void) vulog(LOG_ERR, fmt, args);
	}
	va_end(args);
}


/*
 * Log program errors
 * Calling sequence is
 *	uerror(format, arg1, arg2,...)
 * with zero or more args of types compatible with the associated format
 * specifiers.  For example:
 *
 *         uerror("Inconsistant input %s", input);
 */
#ifdef STDC_ARGS
void
uerror(const char *fmt, ...)
#else
/*VARARGS1*/
void
uerror(fmt, va_alist)
     char *fmt;
     va_dcl
#endif /* STDC_ARGS */
{
    va_list args;
#ifdef STDC_ARGS
	va_start(args ,fmt);
#else
	va_start(args);
#endif /* STDC_ARGS */
	(void) vulog(LOG_ERR, fmt, args);
	va_end(args);
}


/*
 * Log "Normal but significant conditions"
 * Calling sequence is
 *	unotice(format, arg1, arg2,...)
 * with zero or more args of types compatible with the associated format
 * specifiers.  For example:
 *
 *         unotice("Shutting down on signal %s", s_signal(sig));
 */
#ifdef STDC_ARGS
void
unotice(const char *fmt, ...)
#else
/*VARARGS1*/
void
unotice(fmt, va_alist)
     char *fmt;
     va_dcl
#endif /* STDC_ARGS */
{
    va_list args;
#ifdef STDC_ARGS
	va_start(args ,fmt);
#else
	va_start(args);
#endif /* STDC_ARGS */
	(void) vulog(LOG_NOTICE, fmt, args);
	va_end(args);
}


/*
 * Log informational messages
 * Calling sequence is
 *	uinfo(format, arg1, arg2,...)
 * with zero or more args of types compatible with the associated format
 * specifiers.  For example:
 *
 *         uinfo("%s", info->ident);
 */
#ifdef STDC_ARGS
void
uinfo(const char *fmt, ...)
#else
/*VARARGS1*/
void
uinfo(fmt, va_alist)
     char *fmt;
     va_dcl
#endif /* STDC_ARGS */
{
    va_list args;
#ifdef STDC_ARGS
	va_start(args ,fmt);
#else
	va_start(args);
#endif /* STDC_ARGS */
	(void) vulog(LOG_INFO, fmt, args);
	va_end(args);
}


/*
 * Log debugging info
 * Calling sequence is
 *	udebug(format, arg1, arg2,...)
 * with zero or more args of types compatible with the associated format
 * specifiers.  For example:
 *
 *         udebug("entering myproc, arg = %d", arg);
 * This one is a little different than the others in that it behaves diferently
 * when logging to a file than when logging to syslogd.
 */
#ifdef STDC_ARGS
void
udebug(const char *fmt, ...)
#else
/*VARARGS1*/
void
udebug(fmt, va_alist)
     char *fmt;
     va_dcl
#endif /* STDC_ARGS */
{
    va_list args;
	
#ifdef STDC_ARGS
	va_start(args ,fmt);
#else
	va_start(args);
#endif /* STDC_ARGS */

	if(logFilename != NULL && logFd != -1
			&& (LOG_MASK(LOG_PRI(LOG_DEBUG)) & logMask) != 0 )
	{
		char buf[1024];
		/*
		 * When sending debug info to a file (the usual case),
		 * we don't want all the timestamp and pid stuff.
		 * So we format and write it directly from here.
		 */
		size_t cnt;
		buf[0] = '\t'; buf[1] = 0;
		(void) vsprintf(&buf[1], fmt, args);
		cnt = strlen(buf); /* cnt for write() below */
		if(buf[cnt -1] != '\n')
		{
			buf[cnt++] = '\n';
			buf[cnt] = 0;
		}
#if TBUFDIAG
		fputs(buf, stdout);
#endif
		/* output the message to the file, pray all goes well */
		if(write(logFd, buf, cnt) < 0)
			(void) closeulog();
	}
	else
	{
		(void) vulog(LOG_DEBUG, fmt, args);
	}

	va_end(args);
}

void
_uassert(
	const char *ex,
	const char *file,
	int line)
{
	uerror("assertion \"%s\" failed: file \"%s\", line %d\n",
		ex, file, line);
	abort();
}


#ifdef TIRPC
#include <tiuser.h>
extern int t_errno;
extern char *t_errlist[];
extern int t_nerr;

/*
 * Log tli call errors
 * Use where you would want to call t_error(3N).
 * Calling sequence is
 *	terror(format, arg1, arg2,...)
 * with zero or more args of types compatible with the associated format
 * specifiers.
 */
#ifdef STDC_ARGS
void
terror(const char *fmt, ...)
#else
/*VARARGS1*/
void
terror(fmt, va_alist)
     char *fmt;
     va_dcl
#endif /* STDC_ARGS */
{
    va_list args;
#ifdef STDC_ARGS
	va_start(args ,fmt);
#else
	va_start(args);
#endif /* STDC_ARGS */
	if(t_errno != 0)
	{
		char buf[1024];
		int errnum = t_errno;	/* save real t_errno in case we wipe it out */
		char *cp;
		(void) vsprintf(buf, fmt, args);
		for(cp = buf; *cp != 0; cp++) /*EMPTY*/;
		strcat(cp, ": %s");
		ulog(LOG_ERR, buf,
			(errnum > 0 && errnum < t_nerr) ? t_errlist[errnum] : "");
		t_errno = 0;
	}
	else
	{
		(void) vulog(LOG_ERR, fmt, args);
	}
	va_end(args);
}
#endif


/* strip off leading path */
const char *
ubasename(const char *av0)
{
	const char *logident;
#ifdef vms
#define SEP	']'
#endif
#ifndef SEP
#define SEP	'/'
#endif
	if ((logident = strrchr(av0, SEP)) == NULL)
		logident = av0;
	else
	    logident++;
	return logident;
}


/*
 * regerror() in Henry Spencer's regular-expression library calls exit().
 * Because his is not desirable, we override the definition here.
 */
void
regerror(const char *msg)
{
	uerror("regexp(3l): %s", msg);
}


#ifdef  TEST0
#include <stdlib.h>

void
timecheck(void)
{
	time_t now;
	struct tm local[1];
	struct tm gmt[1];
	char cp[32];
	(void) time(&now);
	*local = *(localtime(&now));
	strftime(cp, sizeof(cp), "%b %d %H:%M:%S ", local);
	unotice("Local: %s (%02d)", cp, local->tm_hour);
	*gmt = *(gmtime(&now)); /* may dump core */
	strftime(cp, sizeof(cp), "%b %d %H:%M:%S ", gmt);
	unotice("  UTC: %s (%02d)", cp, gmt->tm_hour);
}

void
ploop(int ii)
{
	ulog(LOG_ALERT, "%d Alert\n", ii);
	ulog(LOG_CRIT, "%d Crit\n", ii);
	ulog(LOG_ERR, "%d Err\n", ii);
	ulog(LOG_WARNING, "%d Warning\n", ii);
	ulog(LOG_NOTICE, "%d Notice\n", ii);
	ulog(LOG_INFO, "%d Info %d \"%m\" %d %s %f %d\n",
		ii, 1, 2, "string", 3.1415, 4 );
	ulog(LOG_DEBUG, "%d Debug\n", ii);
}

void
usage(char *av0)
{
	fprintf(stderr, "Usage: %s [-l logfname]\n", av0);
	exit(1);
}

main(int ac, char *av[])
{
	char *logfname = 0;
	int logfd = -1;

	{
	extern int optind;
	extern int opterr;
	extern char *optarg;
	int ch;

	opterr = 1;

	while ((ch = getopt(ac, av, "l:")) != EOF)
		switch (ch) {
		case 'l':
			logfname = optarg;
			break;
		case '?':
			usage(av[0]);
			break;
		}
	}

	printf("ulog[%d]\n", getpid());

	logfd = openulog("ulog", (LOG_CONS|LOG_PID), LOG_LDM, logfname);
	if(logfd == -1)
		perror("openulog");

	ploop(0);
	(void) setulogmask(LOG_MASK(LOG_NOTICE));
	ploop(1);
	(void) setulogmask(LOG_UPTO(LOG_INFO));
	ploop(2);
	(void) setulogmask(LOG_UPTO(LOG_ALERT));
	ploop(3);
	(void) setulogmask((LOG_UPTO(LOG_DEBUG) & ~(LOG_UPTO(LOG_CRIT))));
	ploop(4);
	{
	int fd = open("/dev/kmem", O_RDWR, 0664);
	if(fd == -1)
		serror("serror %d %d %s %f %d",
			1, 2, "string", 3.1415, 4 );
	}

	timecheck();

	/* cheating */
	logOptions |= LOG_LOCALTIME;

	timecheck();
	closeulog();

	exit(0);
}
#endif

#ifdef TEST1

void
ploop2(ii)
int ii;
{
	unotice("%d run (notice)", ii);
	errno = 13;
	serror("error: filename");
	serror("noerror");
	uerror("error %d %d %s %f %d", 
			1, 2, "string", 3.1415, 4 );
	uinfo("information");
	udebug("debug  %d %d %s %f %d",
		1, 2, "string", 3.1415, 4 );
	unotice("\t%d end (notice)", ii);
}


main()
{
	int logfd = -1;
	logfd = openulog("ulog", (LOG_CONS|LOG_PID), LOG_LOCAL0, 0);
	if(logfd == -1)
		perror("openulog");
	ploop2(1);
	closeulog();
	logfd = openulog("ulog", (LOG_CONS|LOG_PID), LOG_LOCAL0, "logfile");
	(void) setulogmask((LOG_UPTO(LOG_DEBUG) & ~(LOG_MASK(LOG_INFO))));
	if(logfd == -1)
		perror("openulog");
	ploop2(2);
	closeulog();
}
#endif
