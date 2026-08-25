#!/usr/bin/env python3
# 直接解 CCSDS 原始包抽 raw channel(byte0)做直方，跳过时间 solve。
# 复制 blink 的 parse_events + crc_check 语义：
#   每 882B 包: payload=[6:878], 每 8B 一个事件
#   crc_check(row[0..7]) == row[7]&0x0F 才有效
#   type=row[7]&0x30: 0x00/0x20=事件, 0x10=SEC(排除)
#   channel = byte0
#
# flood 判据的离线原型：ch15-19 本该是近空的 gap，flood 小时把它填满，
# 判据落地在 config_guard.rs。需要集群上的 1B 档案。
#
# usage: channel_hist.py <label> <YYYYMMDD> <hour>
import sys, glob
import numpy as np
from astropy.io import fits

CRC_TABLE = np.array([0,3,6,5,12,15,10,9,11,8,13,14,7,4,1,2], dtype=np.uint8)

def find_evt_file(date, hour):
    # 优先 0642(一致的box),否则任意 HE_Evt_Src
    base = "/hxmtfs/data/Archive_tmp/1B/%s/%s" % (date[:4], date)
    for sub in ["0642", "0922", "1686"]:
        fs = glob.glob("%s/%s/*T%02d0000_*.fits" % (base, sub, hour))
        if fs:
            return fs[0]
    # 兜底：扫所有子目录找 HE_Evt_Src
    for f in glob.glob("%s/*/*T%02d0000_*.fits" % (base, hour)):
        try:
            with fits.open(f) as h:
                if any(getattr(hd, "name", "") == "HE_Evt_Src" for hd in h):
                    return f
        except Exception:
            pass
    return None

def hist_file(f):
    with fits.open(f) as h:
        pk = np.asarray(h["HE_Evt_Src"].data["CCSDS"], dtype=np.uint8)  # (N,882)
    if pk.ndim == 1:
        pk = pk.reshape(-1, 882)
    payload = pk[:, 6:878]                       # (N,872)
    ev = payload.reshape(-1, 8)                  # (M,8) 每行一个事件
    b7 = ev[:, 7]
    # 向量化 CRC: 15 个 nibble = byte0hi,byte0lo,...,byte6hi,byte6lo,byte7hi
    crc = np.zeros(ev.shape[0], dtype=np.uint8)
    nibs = []
    for b in range(7):
        nibs.append((ev[:, b] & 0xF0) >> 4)
        nibs.append(ev[:, b] & 0x0F)
    nibs.append((b7 & 0xF0) >> 4)
    for nib in nibs:
        crc = CRC_TABLE[(crc ^ nib)]
    crc_ok = crc == (b7 & 0x0F)
    typ = b7 & 0x30
    is_evt = (typ == 0x00) | (typ == 0x20)       # 事件(非SEC)
    valid = crc_ok & is_evt
    ch = ev[valid, 0]
    return np.bincount(ch, minlength=256)

def main():
    label, date, hour = sys.argv[1], sys.argv[2], int(sys.argv[3])
    f = find_evt_file(date, hour)
    if f is None:
        print("%s: NO-FILE (%s T%02d)" % (label, date, hour)); return
    hc = hist_file(f)
    tot = int(hc.sum())
    lo = int(hc[:20].sum())
    print("%s: total=%d frac0_19=%.3f file=%s" % (label, tot, (lo/tot if tot else 0), f.split("/")[-1]))
    print("ch " + " ".join("%d:%d" % (i, hc[i]) for i in range(46)))

if __name__ == "__main__":
    main()
