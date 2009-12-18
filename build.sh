#!/bin/bash

usage="Usage: ${0##*/} [options] [progs]"
prefix="$HOME/votca"
#mind the spaces
all=" tools csg "
standard=" tools csg "

do_prefix_clean="no"
do_configure="yes"
do_clean="yes"
do_build="yes"
do_install="yes"
do_update="no"

relurl="http://votca.googlecode.com/files/votca-PROG-REL.tar.gz"
rel=""
url="https://PROG.votca.googlecode.com/hg/"
checkout="hg clone"
default="defaut"

extra_conf=""

gromacs="no"

BLUE="[34;01m"
CYAN="[36;01m"
GREEN="[32;01m"
RED="[31;01m"
OFF="[0m"

die () {
  cecho RED "$*" >&2
  exit 1
}

cecho() {
  local opts color=" BLUE CYAN GREEN RED "
  if [ -z "${1##-*}" ]; then
    opts="$1"
    shift
  fi
  [ -z "$2" ] && die "cecho: Missing argumet"
  [ -n "${color//* $1 *}" ] && die "cecho: Unknown color ($color allowed)"
  color=${!1}
  shift
  echo -n ${color}
  echo -ne "$@"
  echo $opts "${OFF}"
}

prefix_clean() {
  cd $prefix || die "Dir: '$prefix not found'"
  files="$(ls -d bin include lib share 2>/dev/null)"
  if [ -z "$files" ]; then 
    echo "Found nothing to clean"
    cd -
    return
  fi
  echo "I will $(cecho RED remove):"
  echo $files
  cecho RED "CTRL-C to stop it"
  countdown 10
  rm -rf $files
  cecho GREEN "Done, hope you are happy now"
  cd -
}

countdown() {
  [ -z "$1" ] && "countdown: Missing argument"
  [ -n "${1//[0-9]}" ] && "countdown: argument should be a number"
  for ((i=$1;i>0;i--)); do
    cecho -n CYAN "$i "
    sleep 1
  done
  echo
}

show_help () {
  cat << eof
This is the votca build utils which builds votca modules
Give multiple programs to build them. Nothing means:$standard
One can build:$all

Please visit: $(cecho BLUE www.votca.org)

The normal sequence of a build is:
- hg clone (if src is not there)
  (use release tarball with --release)
- hg pull + hg update (enable --do-update)
  (stop here with --no-configure)
- bootstrap
- configure
- make clean (disable with --no-clean)
  (stop here with --no-build)    
- make 
- make install (disable with --no-install)

The most recent version can be found at:
$(cecho BLUE http://hg.votca.org/buildutil/raw-file/tip/build.sh)

$usage
OPTIONS:
$(cecho GREEN -h), $(cecho GREEN --help)              Show this help
    $(cecho GREEN --nocolor)           Disable color
    $(cecho GREEN --votca.org)         Use votca.org server instead of googlecode
                        (less reliable)
$(cecho GREEN -d),  $(cecho GREEN --dev)              Use votca developer repository
                        (username/password needed)
    $(cecho GREEN --ccache)            Enable ccache
    $(cecho GREEN --release) $(cecho CYAN REL)       Get Release tarball instead of use hg clone 
$(cecho GREEN -u), $(cecho GREEN --do-update)         Do a update from hg
$(cecho GREEN -c), $(cecho GREEN --clean-out)         Clean out the prefix
    $(cecho GREEN --no-configure)      Stop after update (before bootstrap)
    $(cecho GREEN --conf-opts) $(cecho CYAN OPTS)    Extra configure options
    $(cecho GREEN --no-clean)          Don't run make clean
    $(cecho GREEN --no-build)          Stop before build
    $(cecho GREEN --no-install)        Don't run make install
    $(cecho GREEN --prefix) $(cecho CYAN \<prefix\>)   use prefix
                        Default: $prefix
$(cecho GREEN -g), $(cecho GREEN --gromacs)           Set gromacs stuff base up your \$GMXLDLIB

Examples:  ${0##*/} tools csg
           ${0##*/} -cug --prefix \$PWD/install tools csg
	   ${0##*/} -u
	   ${0##*/} --release 1.0_rc1 tools csg
	
eof
}

# parse arguments

while [ "${1#-}" != "$1" ]; do
 if [ "${1#--}" = "$1" ] && [ -n "${1:2}" ]; then
    #short opt with arguments here: f
    if [ "${1#-[f]}" != "${1}" ]; then
       set -- "${1:0:2}" "${1:2}" "${@:2}"
    else
       set -- "${1:0:2}" "-${1:2}" "${@:2}"
    fi
 fi
 case $1 in
   -h | --help)
    show_help
    exit 0;;
   -c | --clean-out)
    prefix_clean="yes"
    shift 1;;
   -g | --gromacs)
    gromacs="yes"
    shift 1;;
   -u | --do-update)
    do_update="yes"
    shift 1;;
   --no-configure)
   do_configure="no"
    shift 1;;
   --no-clean)
   do_clean="no"
    shift 1;;
   --no-install)
    do_install="no"
    shift 1;;
   --no-build)
    do_build="no"
    shift 1;;
   --prefix)
    prefix="$2"
    shift 2;;
   --conf-opts)
    extra_conf="$2"
    shift 2;;
   --release)
    rel="$2"
    shift 2;;
   --ccache)
    [ -z "$(type ccache)" ] && die "${0##*/}: ccache not found"
    export CXX="ccache ${CXX:=g++}"
    shift;;
   --nocolor)
    unset BLUE CYAN GREEN OFF RED
    shift;;
   --votca.org)
    url="http://hg.votca.org/PROG"
    relurl="http://www.votca.org/downloads/votca-PROG-REL.tar.gz"
    shift;;
   -d | --dev)
    url="http://dev.votca.org/votca/PROG"
    all=" tools csg moo kmc tof testsuite "
    standard=" tools csg moo kmc tof "
    shift 1;;
  *)
   die "Unknown option '$1'"
   exit 1;;
 esac
done

[ -z "$1" ] && set -- $standard
[ -z "$prefix" ] && die "Error: prefix is empty"


CPPFLAGS="-I$prefix/include $CPPFLAGS"
LDFLAGS="-L$prefix/lib $LDFLAGS"

if [ "$gromacs" = "yes" ]; then
  [ -z "$GMXLDLIB" ] && die "Error: GMXLDLIB is not defined"
  LDFLAGS="-L$GMXLDLIB $LDFLAGS"
  CPPFLAGS="-I$GMXLDLIB/../include/gromacs $CPPFLAGS"
fi
export CPPFLAGS LDFLAGS

echo "prefix is '$prefix'"
echo "CPPFLAGS is '$CPPFLAGS'"
echo "LDFLAGS is '$LDFLAGS'"
[ "$prefix_clean" = "yes" ] && prefix_clean

set -e
for prog in "$@"; do
  [ -n "${all//* $prog *}" ] && die "Unknown progamm '$prog', I know$all"

  cecho GREEN "Working on $prog"
  if [ -d "$prog" ] && [ -z "$rel" ]; then
    cecho BLUE "Source dir is already there - skipping checkout"
  elif [ -d "$prog" ] && [ -n "$rel" ]; then
    cecho BLUE "Source dir is already there - skipping download (CTRL-C to stop)"
    countdown 5
  elif [ ! -d "$prog" ] && [ -n "$rel" ]; then
    tmpurl="${relurl//REL/$rel}"
    tmpurl="${tmpurl//PROG/$prog}"
    tarball="${tmpurl##*/}"
    cecho GREEN "Download tarball $tarball from ${tmpurl}"
    [ -f "$tarball" ] && die "Tarball $tarball is already there, remove it first"
    wget "${tmpurl}"
    tardir="$(tar -tzf ${tarball} | sed -e's#/.*$##' | sort -u)"
    [ -z "${tardir//*\n*}" ] && die "Tarball $tarball contains zero or more then one directory, please check by hand"
    [ -e "${tardir}" ] && die "Tarball unpack directory ${tardir} is already there, remove it first"
    tar -xzf "${tarball}"
    mv "${tardir}" "${prog}"
    rm -f "${tarball}"
  else
    cecho BLUE "Doing checkout for $prog from ${url/PROG/$prog} (CTRL-C to stop)"
    countdown 5
    ${checkout} ${url/PROG/$prog} $prog
  fi

  cd $prog
  if [ "$do_update" == "yes" ]; then
    if [ -n "$rel" ]; then
      cecho BLUE "Update of a release tarball doesn't make sense, skipping (CTRL-C to stop)"
      countdown 5
    elif [ -d .hg ]; then
      cecho GREEN "updating hg repository"
      pullpath=$(hg path $default 2> /dev/null || true )
      if [ -z "${pullpath}" ]; then
	pullpath=${url/PROG/$prog}
	cecho BLUE "Could not fetch '$default' pull path, using $pullpath instead (CTRL-C to stop)"
	countdown 5
      else
	cecho GREEN "from $pullpath"
      fi
      hg pull ${pullpath}
      hg update
    else
      cecho BLUE "$prog dir doesn't seem to be a hg repository, skipping update (CTRL-C to stop)"
      countdown 5
    fi
  fi
  cecho GREEN "configuring $prog"
  if [ "$do_configure" == "yes" ]; then
    [ -f bootstrap.sh ] && ./bootstrap.sh 
    ./configure --prefix "$prefix" $extra_conf
  else
    cd ..
    cecho GREEN "done with $prog"
    continue 
  fi
  if [ "$do_clean" == "yes" ]; then
    cecho GREEN "cleaning $prog"
    make clean
  fi
  if [ "$do_build" == "no" ]; then 
    cd .. 
    cecho GREEN "done with $prog"
    continue
  fi
  cecho GREEN "buidling $prog"
  make
  if [ "$do_install" == "yes" ]; then 
    cecho GREEN "installing $prog"
    make install
  fi
  cd ..
  cecho GREEN "done with $prog"
done
set +e

