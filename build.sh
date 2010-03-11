#!/bin/bash
# 
# Copyright 2009 The VOTCA Development Team (http://www.votca.org)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#version 1.0.0 -- 18.12.09 initial version
#version 1.0.1 -- 21.12.09 added --pullpath option
#version 1.0.2 -- 14.01.10 improved clean
#version 1.0.3 -- 20.01.10 better error message in prefix_clean
#version 1.0.4 -- 09.02.10 added --static option 
#version 1.0.5 -- 03.03.10 added pkg-config support

#defaults
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
selfurl="http://votca.googlecode.com/hg/build.sh"
pathname="default"

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
  cecho GREEN "Starting clean out of prefix"
  [ ! -d $prefix ] && cecho BLUE "prefix '$prefix' is not there - skipping" && return 0
  cd $prefix || die "Could change to prefix '$prefix'"
  files="$(ls -d bin include lib share 2>/dev/null)"
  if [ -z "$files" ]; then 
    cecho BLUE "Found nothing to clean"
    cd - > /dev/null
    return
  fi
  echo "I will $(cecho RED remove):"
  echo $files
  cecho RED "CTRL-C to stop it"
  countdown 10
  rm -rf $files
  cecho GREEN "Done, hope you are happy now"
  cd - > /dev/null
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

get_version() {
  sed -ne 's/^#version[[:space:]]*\([^[:space:]]*\)[[:space:]]*-- .*$/\1/p' $1 | sed -n '$p'
}

self_update() {
  [ -z "$(type -p wget)" ] && die "wget is missing"
  new_version="$(wget -qO- "${selfurl}")" || die "self_update: wget fetch from $selfurl failed"
  new_version="$(echo -e "${new_version}" | get_version)"
  [ -z "${new_version}" ] && die "self_update: Could not fetch new version number"
  cecho BLUE "Version of $selfurl is: $new_version"
  old_version="$(get_version $0)"
  cecho BLUE "Local Version: $old_version"
  new_version="${new_version//[^0-9]}"
  old_version="${old_version//[^0-9]}"
  newer=$(awk -v new="$new_version" -v old="$old_version" 'BEGIN{if (new>old){print "yes"}}')
  if [ "$newer" = "yes" ]; then
    cecho RED "I will try replace myself now with $selfurl (CTRL-C to stop)"
    countdown 5
    wget -O "${0}" "${selfurl}"
  else
    cecho GREEN "No updated needed"
  fi
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
- bootstrap (if needed)
- configure
- make clean (disable with --no-clean)
  (stop here with --no-build)    
- make 
- make install (disable with --no-install)

The most recent version can be found at:
$(cecho BLUE $selfurl)

$usage
OPTIONS (last overwrites previous):
$(cecho GREEN -h), $(cecho GREEN --help)              Show this help
$(cecho GREEN -v), $(cecho GREEN --version)           Show version
    $(cecho GREEN --debug)             Enable debug mode
    $(cecho GREEN --nocolor)           Disable color
    $(cecho GREEN --votca.org)         Use votca.org server instead of googlecode
                        (less reliable)
    $(cecho GREEN --selfupdate)        Do a self update (EXPERIMENTAL)
$(cecho GREEN -d), $(cecho GREEN --dev)               Use votca developer repository
                        (username/password needed)
    $(cecho GREEN --ccache)            Enable ccache
    $(cecho GREEN --static)            Build static executables
    $(cecho GREEN --release) $(cecho CYAN REL)       Get Release tarball instead of using hg clone 
$(cecho GREEN -u), $(cecho GREEN --do-update)         Do a update of the sources from pullpath $pathname
                        or the votca server as fail back
    $(cecho GREEN --pullpath) $(cecho CYAN NAME)     Changes the name of the path to pull from
                        Default: $pathname (Also see 'hg paths --help')
$(cecho GREEN -c), $(cecho GREEN --clean-out)         Clean out the prefix (DANGEROUS)
    $(cecho GREEN --no-configure)      Stop after update (before bootstrap)
    $(cecho GREEN --conf-opts) $(cecho CYAN OPTS)    Extra configure options
    $(cecho GREEN --no-clean)          Don't run make clean
    $(cecho GREEN --no-build)          Stop before build
    $(cecho GREEN --no-install)        Don't run make install
    $(cecho GREEN --prefix) $(cecho CYAN \<prefix\>)   use prefix
                        Default: $prefix
$(cecho GREEN -g), $(cecho GREEN --gromacs)           Set gromacs stuff base up your \$GMXLDLIB

Examples:  ${0##*/} tools csg
           ${0##*/} -dcug --prefix \$PWD/install tools csg
	   ${0##*/} -u
	   ${0##*/} --release 1.0_rc1 tools csg
	   ${0##*/} --dev --help
	
eof
}

# parse arguments

shopt -s extglob
while [ "${1#-}" != "$1" ]; do
 if [ "${1#--}" = "$1" ] && [ -n "${1:2}" ]; then
    #short opt with arguments here: f
    if [ "${1#-[r]}" != "${1}" ]; then
       set -- "${1:0:2}" "${1:2}" "${@:2}"
    else
       set -- "${1:0:2}" "-${1:2}" "${@:2}"
    fi
 fi
 case $1 in
   --debug)
    set -x
    shift ;;
   -h | --help)
    show_help
    exit 0;;
   -v | --version)
    echo "${0##*/}, version $(get_version $0)"
    exit 0;;
   --selfupdate)
    self_update
    exit $?;;
   -c | --clean-out)
    prefix_clean="yes"
    shift 1;;
   -g | --gromacs)
    gromacs="yes"
    shift 1;;
   -u | --do-update)
    do_update="yes"
    shift 1;;
   --pullpath)
    pathname="$2"
    shift 2;;
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
    extra_conf="${extra_conf}$2 "
    shift 2;;
   --static)
    extra_conf="${extra_conf}--enable-all-static "
    shift ;;
   --release)
    rel="$2"
    [ -z "${rel//[1-9].[0-9]?(_rc[1-9]?([0-9]))}" ] || \
      die "--release option needs an argument of the form X.X{_rcXX}"
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
    selfurl="http://hg.votca.org/buildutil/raw-file/tip/build.sh"
    shift;;
   -d | --dev)
    url="http://dev.votca.org/votca/PROG"
    #wget would needs auth
    selfurl="http://dev.votca.org/votca/buildutil/raw-file/tip/build.sh"
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

export PKG_CONFIG_PATH="$prefix/lib/pkgconfig${PKG_CONFIG_PATH:+:}${PKG_CONFIG_PATH}"
if [ "$gromacs" = "yes" ]; then
  [ -z "$GMXLDLIB" ] && die "Error: GMXLDLIB is not defined"
  if [ -f "$GMXLDLIB/pkgconfig/libgmx.pc" ]; then
    export PKG_CONFIG_PATH="$GMXLDLIB/pkgconfig:${PKG_CONFIG_PATH}"
  else
    export GMX_LIBS="-L$GMXLDLIB -lgmx"
    export GMX_CFLAGS="-I$GMXLDLIB/../include/gromacs"
  fi
fi

echo "prefix is '$prefix'"
[ -n "$CPPFLAGS" ] && echo "CPPFLAGS is '$CPPFLAGS'"
[ -n "$LDFLAGS" ] && echo "LDFLAGS is '$LDFLAGS'"
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
    [ -z "$(type -p wget)" ] && die "wget is missing"
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
    hg clone ${url/PROG/$prog} $prog
  fi

  cd $prog
  if [ "$do_update" == "yes" ]; then
    if [ -n "$rel" ]; then
      cecho BLUE "Update of a release tarball doesn't make sense, skipping (CTRL-C to stop)"
      countdown 5
    elif [ -d .hg ]; then
      cecho GREEN "updating hg repository"
      pullpath=$(hg path $pathname 2> /dev/null || true )
      if [ -z "${pullpath}" ]; then
	pullpath=${url/PROG/$prog}
	cecho BLUE "Could not fetch pull path '$pathname', using $pullpath instead (CTRL-C to stop)"
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
  if [ "$do_configure" == "yes" ]; then
    if [ -f bootstrap.sh ]; then 
      cecho GREEN "bootstraping $prog"
      ./bootstrap.sh
    fi
    cecho GREEN "configuring $prog"
    echo configure --prefix "$prefix" $extra_conf
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

