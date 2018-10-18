#include <stdio.h>
#include <stdlib.h> 
#include <errno.h>
#include <string.h>
#include <netinet/in.h>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <arpa/inet.h>
#include <pwd.h>
#include <syslog.h>

#include "mond.h"

/* For now, the only way to turn on verbosity is at compile time.
   verbose mode generates syslog messages with priority LOG_DEBUG.
   These may or may not show up somewhere (/var/log/messages maybe?)
   depending on the settings in /etc/syslog.conf.  */
int verbose = 1;

/* safely copy the 'new' string to 'to', without writing on or past 'end'.
   If the initial value of 'to' or 'new' is NULL or if 'to' is greater
   than end, do nothing and return NULL.  Otherwise, always ensure
   that the result is NUL-terminated, and return a pointer to the
   trailing NUL (note the return is very different from the standard
   strcpy's).  Note that sufficiently bad input can still segfault (e.g.,
   if new is 0x1), but at least it's impossible to write outside the
   'sandbox' between to and end.
*/
char *strcpy_safe(char *to, const char *new, const char *end){
  char c;

  if( to==NULL || new==NULL || to>=end ) return NULL;
  while( to<end && (c=*new++) ) *to++ = c;
  if( to == end ) --to;
  *to = '\0';
  return to;
}

int
main(int argc, char **argv)
{
  int sockfd;			/* listen on sock_fd */
  struct sockaddr_in my_addr;   /* my address information */
  struct sockaddr_in their_addr; /* connector's address information */
  int sin_size;
  int pid;
  char buffer[MAXBUFFER], input[MAXBUFFER];
  char *end;
  FILE *fp;

  /* Immediately return control to the caller. */
  pid = fork();
  if (pid != 0) exit(0);

  /* Don't disturb stdin, stdout, stderr */
  close(0);
  close(1);
  close(2);
  openlog("mond", 0, LOG_DAEMON);
  if( geteuid()==0 ){		/* don't run as root. */
    struct passwd *pw;
    if( (pw=getpwnam("nobody")) == NULL ){
      syslog(LOG_ERR, "can't find pw entry for nobody: %m");
      exit(2);
    }
    if( setuid(pw->pw_uid) < 0 ){
      syslog(LOG_ERR, "Could not setuid to %d: %m", pw->pw_uid);
      exit(3);
    }
  }

  if ((sockfd = socket(AF_INET, SOCK_DGRAM, 0)) == -1) {
    syslog(LOG_ERR, "socket: %m");
    exit(1);
  }

  memset(&my_addr, 0, sizeof(my_addr));    /* zero everything */
  my_addr.sin_family = AF_INET;     /* host byte order */
  my_addr.sin_port = htons(PORT); /* short, network byte order */
  my_addr.sin_addr.s_addr = INADDR_ANY;

  if (bind(sockfd, (struct sockaddr *)&my_addr, sizeof(my_addr)) == -1) {
    syslog(LOG_ERR, "bind: %m");
    exit(1);
  }
  end = buffer + MAXBUFFER;
  while(1) {
    int len;
    char *outp;
    sin_size = sizeof(their_addr);
    if( (len=recvfrom(sockfd, input, MAXBUFFER, 0, &their_addr, &sin_size)) <= 0 ){
      syslog(LOG_ERR, "recvfrom: %m");
      continue;
    }
    buffer[0] = '\0';		/* you can't be too careful */
    outp = strcpy_safe(buffer, input, end);

    if( verbose )
      syslog(LOG_DEBUG, "got connection from %s (port=%d).  initial string: %s", inet_ntoa(their_addr.sin_addr), htons(their_addr.sin_port), buffer);

    outp = strcpy_safe(outp, " ", end);
    if ((fp=fopen("/proc/loadavg","r")) == NULL){
      syslog(LOG_ERR, "fopen(/proc/loadavg): %m");
      continue;			/* exit(1)? */
    }

    /* read /proc/loadavg */
    if( fgets(input,MAXBUFFER,fp) == NULL ){
      syslog(LOG_ERR, "could not read from /proc/loadavg: %m");
      continue;
    }
    outp = strcpy_safe(outp, input, end);
    if( outp && outp[-1] == '\n' ) outp[-1] = ' ';
    fclose(fp);

    /* read /proc/meminfo */
    if ((fp=fopen("/proc/meminfo","r")) == NULL){
      syslog(LOG_ERR, "fopen(/proc/meminfo): %m");
      continue;			/* exit(1)? */
    }

    fgets(input,MAXBUFFER,fp); /* discard the first line */

    if( fgets(input,MAXBUFFER,fp) == NULL ){/* read the second line */
      syslog(LOG_ERR, "could not read from /proc/meminfo: %m");
      continue;
    }
    if( strlen(input) < 5 ){
      syslog(LOG_ERR, "short read from line 2 of /proc/meminfo");
      continue;
    }
    outp = strcpy_safe(outp, input+4, end); /* skip 'Mem:' */
    if( outp && outp[-1] == '\n' ) outp[-1] = ' ';

    if( fgets(input,MAXBUFFER,fp) == NULL ){/* read the third line */
      syslog(LOG_ERR, "could not read from /proc/meminfo: %m");
      continue;
    }
    if( strlen(input) < 6 ){
      syslog(LOG_ERR, "short read from line 3 of /proc/meminfo");
      continue;
    }
    outp = strcpy_safe(outp, input+5, end); /* skip 'Swap:' */
    if( outp && outp[-1] == '\n' ) outp[-1] = ' ';
    fclose(fp);
		
    /* read /proc/meminfo */
    if ((fp=fopen("/proc/sys/dev/sensors/w83782d-i2c-0-2e/temp3","r")) == NULL){
      syslog(LOG_ERR, "fopen(/proc/sys/dev/sensors/w83782d-i2c-0-2e/temp3): %m");
      continue;			/* exit(1)? */
    }
    if( fgets(input,MAXBUFFER,fp) == NULL ){
      syslog(LOG_ERR, "could not read from /proc/meminfo: %m");
      continue;
    }
    outp = strcpy_safe(outp, input+10, end);
    fclose(fp);

    if( outp == NULL || outp >= end ){
      syslog(LOG_ERR, "buffer overflow detected, outp=%p, end=%p",outp,end);
      continue;
    }

    len = (outp - buffer) + 1;
    if( verbose )
      syslog(LOG_DEBUG, "send %d bytes: %s", len, buffer);
    if( sendto(sockfd, buffer, len, 0, &their_addr, sin_size) != len ){
      syslog(LOG_ERR, "sendto: %m");
    }
  }
  /* NOTREACHED */
}
