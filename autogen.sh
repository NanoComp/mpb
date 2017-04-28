#!/bin/sh

configure_args=""

while test $# -ge 1; do
    case $1 in
	--prefix=*) configure_args="$configure_args $1" ;;
	--enable-*) configure_args="$configure_args $1" ;;
	--disable-*) configure_args="$configure_args $1" ;;
	--with-*) configure_args="$configure_args $1" ;;
	--without-*) configure_args="$configure_args $1" ;;
	*) echo "unknown argument $1"; exit 1 ;;
    esac
    shift
done

# paranoia: sometimes autoconf doesn't get things right the first time
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force

./configure --enable-maintainer-mode $configure_args
