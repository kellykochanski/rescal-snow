#!/usr/bin/env bash
#
# Script to convert matplotlib pngs to animated gif
# Arguments: gif name followed by a list of PNG images to turn into a gif
#
# Author: Izaak Beekman <izaak@izaakbeekman.com>
# Date: 2019-10-14
#
# The author waives all rights to the source code in this script and contributes
# it to the rescale-snow project in accordance with the projects GPLv3 license.
# Please see LICENSE and NOTICE in the project root directory for further information.
#
# Obligatory CYA:
#
# THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
# APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
# HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT
# WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND
# PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE
# DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR
# CORRECTION.
#
# IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN
# WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES
# AND/OR CONVEYS THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR
# DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL
# DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM
# (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING RENDERED
# INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE
# OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN IF SUCH
# HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGES.

set -o errexit
set -o nounset
set -o pipefail

if (( $# < 2 )) ; then
    echo ""
    echo "Usage:"
    echo "    $0 <output_name.gif> <png_img1> [<png_img2> [ ... <png_imgN>]]"
    echo "    Requires ImageMagick to be installed."
    echo ""
    exit 1
fi

if ! command -v convert > /dev/null 2>&1 ; then
    echo "You must have imagemagick's \`convert\` command installed and on your \$PATH!"
fi

OUTFILE="$1"
shift

echo "GIF will be written to $OUTFILE, from inputs:"
printf "%s\n" "${@}"

# Create a temporary directory to hold cleaned PNGs
tmp_dir=$(mktemp -d -t rescal-images-XXXXXXXXXX)
# echo "Temporary work directory: ${tmp_dir}"
# ls -ld "${tmp_dir}"

# Cleanup temp directory on exit or abort/user interrupt
function cleanup {
    rm -rf "${tmp_dir}" || true
}
trap cleanup EXIT

# Normalize png images, & create array of cleaned up pngs
declare -a IMG_TO_PROCESS
for img in "$@" ; do
    img_no_path="$(basename "$img")"
    convert -flatten -background White "${img}" "${tmp_dir}/${img_no_path}"
    IMG_TO_PROCESS+=("${tmp_dir}/${img_no_path}")
done

#ls "${tmp_dir}"

# Create the gif from the cleaned up PNGs
convert -delay 20 "${IMG_TO_PROCESS[@]}" -background White -loop 0 "${OUTFILE}"
