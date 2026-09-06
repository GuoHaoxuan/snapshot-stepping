#!/bin/bash
# 在计算节点上冷编译（登录节点负载 30–40、重进程会被掐；/workfs2 在计算节点只读，所以源码放 scratchfs2）。
export PATH=$HOME/.cargo/bin:$PATH
export CARGO_HOME=/scratchfs2/gecam/guohx/.cargo RUSTUP_HOME=/scratchfs2/gecam/guohx/.rustup
cd /scratchfs2/gecam/guohx/build_src && cargo build --release --offline -p blink > /scratchfs2/gecam/guohx/gridrun/farm_logs/build.out 2>&1
code=$?; echo "BUILD_EXIT=$code" >> /scratchfs2/gecam/guohx/gridrun/farm_logs/build.out
[ $code -eq 0 ] && cp /scratchfs2/gecam/guohx/build_src/target/release/blink /scratchfs2/gecam/guohx/gridrun/blink && echo "binary copied" >> /scratchfs2/gecam/guohx/gridrun/farm_logs/build.out
