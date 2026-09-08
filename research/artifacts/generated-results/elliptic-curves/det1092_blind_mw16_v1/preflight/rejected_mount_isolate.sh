#!/bin/bash
# Invoked only inside unshare -Urnmipf --mount-proc. No repository mounted.
set -euo pipefail
root=$1
fixture=$2
worker=$3
output=$4
mode=$5
runtime=/home/royvanrijn/.local/share/jacobian-sage-10.9
mount --make-rprivate /
mkdir -p "$root"/{usr,lib,lib64,bin,etc,dev,proc,tmp,input,work,output} "$root$runtime"
for p in /usr /lib /lib64 /bin "$runtime"; do
    mount --bind "$p" "$root$p"
    mount -o remount,bind,ro "$root$p"
done
for p in /dev/null /dev/urandom /dev/random; do
    touch "$root$p"
    mount --bind "$p" "$root$p"
done
mount -t proc proc "$root/proc"
printf 'root:x:0:0:root:/tmp:/bin/bash\n' > "$root/etc/passwd"
printf 'root:x:0:\n' > "$root/etc/group"
touch "$root/work/worker.py"
mount --bind "$worker" "$root/work/worker.py"
mount -o remount,bind,ro "$root/work/worker.py"
if [ "$mode" != selftest ]; then
    touch "$root/input/fixture.json"
    mount --bind "$fixture" "$root/input/fixture.json"
    mount -o remount,bind,ro "$root/input/fixture.json"
fi
mount --bind "$output" "$root/output"
cd "$root"
args=()
if [ "$mode" = selftest ]; then args=(--selftest); fi
exec chroot "$root" /usr/bin/env -i HOME=/tmp PATH="$runtime/bin:/usr/bin:/bin" PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 "$runtime/bin/python" /work/worker.py "${args[@]}"
