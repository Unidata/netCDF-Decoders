/*
 *   Copyright 1993, University Corporation for Atmospheric Research
 *   See ../COPYRIGHT file for copying and redistribution conditions.
 */
/* $Id: inetutil.c,v 1.3 2003/02/27 21:47:47 rkambic Exp $ */

/* 
 * Miscellaneous functions to make dealing with internet addresses easier.
 */

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <sys/param.h>
#include <string.h>
#include <stdlib.h>
#include "inetutil.h"
#include "ulog.h"

/* error number for gethostby...() functions */
#ifndef h_errno
extern int	h_errno;	
#endif

#ifndef MAXHOSTNAMELEN
#define	MAXHOSTNAMELEN	256
#endif /* ! MAXHOSTNAMELEN */

/*
 * Return a string indicating the problem with one of the gethostby...()
 * functions.
 */
static const char*
host_err_str(void)
{
    static char	msgstr[200];

    switch (h_errno)
    {
    case 0:
	msgstr[0] = 0;
	break;
#ifdef HOST_NOT_FOUND
    case HOST_NOT_FOUND:
	(void) strcpy(msgstr, "no such host is known");
	break;
#endif
#ifdef TRY_AGAIN
    case TRY_AGAIN:
	(void) strcpy(msgstr,
	    "local server did not receive authoritative response");
	break;
#endif
#ifdef NO_RECOVERY
    case NO_RECOVERY:
	(void) strcpy(msgstr, "nonrecoverable error");
	break;
#endif
#ifdef NO_ADDRESS
    case NO_ADDRESS:
#endif
	(void) strcpy(msgstr, "valid name has no IP address");
	break;
    default:
	sprintf(msgstr, "h_errno = %d", h_errno);
    }

    return msgstr;
}


/*
 * convenience wrapper around gethostname(2)
 * !NO_INET_FQ_KLUDGE ==> try to make sure it is "fully qualified"
 */
char *
ghostname(void)
{

	static char hostname[MAXHOSTNAMELEN+1];

	if (hostname[0])
		return hostname;

	/*
	 * Since the ldm programs require fully qualified
	 * hostnames in an internet environment AND users
	 * often don't have control over the system admin
	 * conventions, we allow override of the 
	 * hostname with an environment variable.
	 * This is meant to be the fully qualified
	 * hostname of the current host.
	 */
	{
		char *cp;
		cp = getenv("LDMHOSTNAME");
		if(cp != NULL)
		{
			strncpy(hostname, cp, MAXHOSTNAMELEN);
			return hostname;
		}
	}

	if(gethostname(hostname, MAXHOSTNAMELEN) < 0)
		return NULL;
#ifndef NO_INET_FQ_KLUDGE
	if(strchr(hostname, '.') == NULL)
	{
		/* gethostname return not fully qualified */
		struct hostent *hp;
		hp = gethostbyname(hostname);
		if(hp != NULL && hp->h_addrtype == AF_INET) 
		{
			/* hopefully hp->h_name is fully qualified */
			strncpy(hostname, hp->h_name, MAXHOSTNAMELEN);
		}
	}
	/* 
	 * On some systems, neither gethostname() nor
	 * the hp->h_name is fully qualified.
	 * If you can't shoot the Systems Administrator and fix it,
	 * hardwire the trailing path here.
	 * (Uncomment and replace ".unversity.edu" with your domain.)
	 */
/* #define HARDWIRED_LOCAL_DOMAIN ".unversity.edu" */
#ifdef HARDWIRED_LOCAL_DOMAIN
	if(strchr(hostname, '.') == NULL)
	{
		strcat(hostname, HARDWIRED_LOCAL_DOMAIN);
	}
#endif /* HARDWIRED_LOCAL_DOMAIN */
#endif
	return hostname;
}


/*
 * Return a string identifying the internet host referred to
 * by paddr. If the hostname lookup fails, returns "dotted quad"
 * form of the address.
 */
char *
hostbyaddr(
	struct sockaddr_in *paddr
)
{
	struct hostent *hp;
	char *otherguy = NULL;
	struct in_addr addr;

	addr.s_addr = paddr->sin_addr.s_addr;

	hp = gethostbyaddr((char *) &addr,
		sizeof (addr), AF_INET);
	if(hp != NULL)
		otherguy = hp->h_name;
	else
		otherguy = inet_ntoa(paddr->sin_addr);
	return otherguy;
}


/*
 * get the sockaddr_in corresponding to 'hostname' (name or NNN.NNN.NNN.NNN)
 */
int
addrbyhost(
	const char *hostname,
	struct sockaddr_in *paddr /* modified on return */
)
{
	
	paddr->sin_addr.s_addr = inet_addr(hostname); /* handle number form */
	if((int) (paddr->sin_addr.s_addr) == -1)
	{
		struct hostent *hp;
		hp = gethostbyname(hostname);
		if(hp == NULL || hp->h_addrtype != AF_INET) 
		{
			return -1;
		}
		(void) memcpy((char *)&paddr->sin_addr, hp->h_addr, (size_t)hp->h_length);
	}
	paddr->sin_family= AF_INET;
	paddr->sin_port= 0;
	(void) memset(paddr->sin_zero, 0, sizeof (paddr->sin_zero));
	return 0;
}


char *
s_sockaddr_in(
	struct sockaddr_in *paddr
)
{
	static char buf[64];
	(void) sprintf(buf,
		"sin_port %5d, sin_addr %s",
		paddr->sin_port,
		inet_ntoa(paddr->sin_addr));
	return buf;
}


/*
 * Puts the address of the current host into *paddr
 * Returns 0 on success, -1 failure
 */
int
gethostaddr_in(
	struct sockaddr_in *paddr
)
{
	char hostname[MAXHOSTNAMELEN];
	struct hostent *hp;

	if(gethostname(hostname,MAXHOSTNAMELEN) == -1)
		return -1;

	hp = gethostbyname(hostname);
	if(hp == NULL || hp->h_addrtype != AF_INET) 
		return -1;
	
	(void) memcpy(&paddr->sin_addr, hp->h_addr, (size_t)hp->h_length);
	paddr->sin_family= AF_INET;
	paddr->sin_port= 0;
	(void) memset(paddr->sin_zero, 0, sizeof (paddr->sin_zero));

	return 0;
}


/*
 * Return the well known port for (servicename, proto)
 * or -1 on failure.
 */
int
getservport(
	const char *servicename,
	const char *proto
)
{
	struct servent *se;
	se = getservbyname(servicename, proto);
	if(se == NULL)
		return -1;
	/* else */
	return se->s_port;
}


/*
 * Attempt to connect to a unix domain socket.
 * Create & connect.
 * Returns (socket) descriptor or -1 on error.
 */
int
usopen(
	const char *name /* name of socket */
)
{
	int sock = -1;
	struct sockaddr addr;	/* AF_UNIX address */

	sock = socket(AF_UNIX, SOCK_DGRAM, 0);
	if(sock == -1)
		return -1;
	/* else */

	addr.sa_family = AF_UNIX;
	(void) strncpy(addr.sa_data, name, sizeof(addr.sa_data));

	if (connect(sock, &addr, sizeof(addr)) == -1)
	{
		(void) close(sock);
		return -1;
	}
	/* else */

	return sock;
}


/*
 * Attempt to connect to a internet domain udp socket.
 * Create & connect.
 * Returns (socket) descriptor or -1 on error.
 */
int
udpopen(
	const char *hostname, 
	const char *servicename
)
{
	int sock = -1;
	int port = -1;
	struct sockaddr_in addr;	/* AF_INET address */

	sock = socket(AF_INET, SOCK_DGRAM, 0);
	if(sock == -1)
		return -1;
	/* else */

	if(addrbyhost(hostname, &addr) == -1)
	{
		(void) close(sock);
		return -1;
	}
	/* else */

	if((port = getservport(servicename, "udp")) == -1)
	{
		(void) close(sock);
		return -1;
	}
	/* else */
	addr.sin_port = (unsigned short) port;

	if (connect(sock, (struct sockaddr *)&addr, sizeof(addr)) == -1)
	{
		(void) close(sock);
		return -1;
	}
	/* else */

	return sock;
}


/*
 * Macro for rounding-up the positive value x to the nearest multiple of n:
 */
#undef	ROUNDUP
#define	ROUNDUP(x,n)	((x % n) ? (x + (n - (x % n))) : x)


/*
 * Return a new (allocated) host entry.
 */
static
struct hostent*
hostent_new(
    const char		*name
)
{
    struct hostent	*new = NULL;
    struct hostent	*entry;

    /*
     * Retrieve the host's entry.
     */
    entry = gethostbyname(name);
    if (NULL == entry)
	uerror("Couldn't get information on host %s: %s", name, host_err_str());
    else
    {
	int		num_aliases;
	int		num_addr;
	char		**from_alias;
	char		**to_alias;
	char		*cp;
	size_t		nbytes;
	size_t		h_name_off;
	size_t		h_aliases_off;
	size_t		addrs_off;
	size_t		h_addr_list_off;
	struct in_addr	**from_addr;
	struct in_addr	**to_addr;
	struct in_addr	*addrp;

	/*
	 * Compute the size requirements and offsets for the new host entry.
	 */

	/* Basic size of host entry structure: */
	nbytes = sizeof(struct hostent);

	/* Offset and size of official name string: */
	h_name_off = nbytes;
	nbytes += strlen(entry->h_name) + 1;

	/* Offset and length of aliases: */
	nbytes = ROUNDUP(nbytes, sizeof(char*));
	h_aliases_off = nbytes;
	for (from_alias = entry->h_aliases; NULL != *from_alias; from_alias++)
	     nbytes += strlen(*from_alias) + 1;
	num_aliases = from_alias - entry->h_aliases;
	nbytes += sizeof(char*) * (num_aliases + 1);

	/* Offset and length of addresses: */
	nbytes = ROUNDUP(nbytes, sizeof(struct in_addr*));
	h_addr_list_off = nbytes;
	for (from_addr = (struct in_addr**)entry->h_addr_list;
			NULL != *from_addr;
			from_addr++)
	    /* EMPTY */;
	num_addr = from_addr - (struct in_addr**)entry->h_addr_list;
	nbytes += sizeof(struct in_addr*) * (num_addr + 1);
	nbytes = ROUNDUP(nbytes, sizeof(struct in_addr));
	addrs_off = nbytes;
	nbytes += sizeof(struct in_addr) * num_addr;

	/*
	 * Allocate the new host entry.
	 */
	new = (struct hostent *) malloc(nbytes);
	if (NULL == new)
	    serror(
	    "Couldn't allocate %lu bytes for information on host \"%s\"", 
		   (unsigned long) nbytes, name);
	else
	{
	    /* Copy non-pointer members. */
	    new->h_addrtype = entry->h_addrtype;
	    new->h_length = entry->h_length;

	    /* Copy official host name. */
	    new->h_name = (char *) new + h_name_off;
	    (void) strcpy(new->h_name, entry->h_name);

	    /* Copy aliases. */
	    new->h_aliases = (char**)((char*)new + h_aliases_off);
	    cp = (char *) (new->h_aliases + num_aliases);
	    for (from_alias = entry->h_aliases, to_alias = new->h_aliases;
		 NULL != *from_alias;
		 from_alias++, to_alias++)
	    {
		*to_alias = cp;
		(void) strcpy(*to_alias, *from_alias);
		cp += strlen(*to_alias) + 1;
	    }
	    *to_alias = NULL;

	    /* Copy addresses. */
	    new->h_addr_list = (char**)((char*)new + h_addr_list_off);
	    addrp = (struct in_addr*)((char*)new + addrs_off);
	    for (from_addr = (struct in_addr**)entry->h_addr_list,
		    to_addr = (struct in_addr**)new->h_addr_list;
		 NULL != *from_addr;
		 from_addr++, to_addr++)
	    {
		*to_addr = addrp++;
		**to_addr = **from_addr;
	    }
	    *to_addr = NULL;
	}					/* new host entry allocated */
    }						/* host entry retrieved */

    return new;
}


/*
 * Compare two (possibly fully-qualified) hostnames.  Indicate if they
 * refer to the same host.  If one of them isn't fully-qualified, then
 * assume it's in the same domain as the other.
 *
 * Returns:
 *	0	Different hosts
 *	1	Same host
 */
static int
same_host(
    const char	*name1,
    const char	*name2
)
{
    return (name1 == name2) ||
	   (strcmp(name1, name2) == 0) ||
	   (strstr(name1, name2) == name1 && name1[strlen(name2)] == '.') ||
	   (strstr(name2, name1) == name2 && name2[strlen(name1)] == '.');
}


/*
 * Attempt to determine if "remote" is the same as this host.
 * Could be lots smarter...
 */
int
isMe(
	const char *remote
)
{
	static char *names[] = {
		"localhost",
		"loopback",
		NULL /* necessary terminator */
	};
	char *me;
	char **npp;
	static struct hostent *hp;

	/* Check `local host' aliases. */
	for (npp = names; *npp != NULL; npp++)
		if (same_host(remote, *npp))
			return !0;

	me = ghostname();
	if (me == NULL)
		return 0;

	/* Check my nominal hostname. */
	if (same_host(me, remote))
		return 1;

	/* Cache host information on myself. */
	if (NULL == hp)
		hp = hostent_new(me);

	/* Check my aliases. */
	if (hp != NULL)
	{
		for(npp = hp->h_aliases; *npp != NULL; npp++)
			if (same_host(*npp, remote))
				return 1;
	}

	return 0;
}


#ifndef TIRPC
/*
 * Create a socket of type "udp" or "tcp" and bind it
 * to port.
 * Return the socket or -1 on error.
 */
int
sockbind(
	const char *type,
	unsigned short port
)
{
	int sock = -1;
	struct sockaddr_in addr;
	size_t len = sizeof(struct sockaddr_in);

	if(type == NULL)
		return -1;

	if(*type == 't')
		sock = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
	else if(*type == 'u')
		sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);

	if(sock == -1) 
		return -1;

	/*
	 * Eliminate problem with EADDRINUSE for reserved socket.
	 * We get this if an upstream data source hasn't tried to
	 * write on the other and we are in FIN_WAIT_2
	 */
	if(*type == 't')
	{
		int on = 1;
		(void) setsockopt(sock, SOL_SOCKET, SO_REUSEADDR,
			(char *)&on, sizeof(on));
	}

	(void) memset((char *)&addr, 0, len);
	addr.sin_family = AF_INET;
	addr.sin_addr.s_addr = INADDR_ANY;
	addr.sin_port = htons(port);

	if(bind(sock, (struct sockaddr *)&addr, len) < 0)
	{
		(void) close(sock);		
		return -1;
	}

	return sock;
}
#else
/*
 * TLI version of above function.
 * Create a TLI transport endpoint of type "udp" or "tcp"
 * and bind it to port.
 * Return the descriptor or -1 on error.
 * Only tested on SunOS 5.x
 */
#include <tiuser.h>
#include <fcntl.h>

int
sockbind(
	const char *type,
	unsigned short port
)
{
	int sock = -1;
	struct sockaddr_in sin_req;
	struct t_bind *req, *ret;
	extern void terror();

	if(type == NULL)
		return -1;

	if(*type == 't')
		sock = t_open("/dev/tcp", O_RDWR, NULL);
	else if(*type == 'u')
		sock = t_open("/dev/udp", O_RDWR, NULL);

	if((sock == -1) || ( t_getstate(sock) != T_UNBND) )
	{
		terror("sockbind: t_open");
		goto err0;
	}

	req = (struct t_bind *)t_alloc(sock, T_BIND, T_ADDR);
	if(req == NULL)
	{
		terror("sockbind: t_alloc req");
		goto err1;
	}
	ret = (struct t_bind *)t_alloc(sock, T_BIND, T_ADDR);
	if(ret == NULL)
	{
		terror("sockbind: t_alloc ret");
		goto err2;
	}

	(void) memset((char *)&sin_req, 0, sizeof(sin_req));
	sin_req.sin_family = AF_INET;
	sin_req.sin_addr.s_addr = INADDR_ANY;
	sin_req.sin_port = htons(port);

	(void) memcpy(req->addr.buf, (char *)&sin_req, sizeof(sin_req));
	req->addr.len = sizeof(sin_req);
	req->qlen = 32; /* rpc_soc.c uses 8 */

	if(t_bind(sock, req, ret) < 0)
	{
		terror("sockbind: t_bind");
		goto err3;
	}
	if(memcmp(req->addr.buf, ret->addr.buf, ret->addr.len) != 0)
	{
		uerror("sockbind: memcmp: t_bind changed address");
	}

	(void) t_free((char *)req, T_BIND);
	(void) t_free((char *)ret, T_BIND);
	return sock;

err3 :
	(void) t_free((char *)ret, T_BIND);
err2 :
	(void) t_free((char *)req, T_BIND);
err1 :
	(void) close(sock);		
err0 :
	return -1;
}

#endif /* !TIRPC */
