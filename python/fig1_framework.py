#!/usr/bin/env python3
"""Generate the study framework figure as a vector PDF.

This figure contains no empirical estimates. Layout parameters are controlled by
config/figures.yml so the visual can be generated consistently without manual editing.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import yaml
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle, PathPatch
from matplotlib.path import Path as MplPath


def brace(ax, x, y0, y1, side="right", lw=1.3, color="black"):
    """Small vertical curly brace in axes coordinates."""
    s = 1 if side == "right" else -1
    dx = 0.012 * s
    ym = (y0 + y1) / 2
    verts = [
        (x, y1), (x+dx, y1), (x+dx, ym+0.02), (x+2*dx, ym),
        (x+dx, ym-0.02), (x+dx, y0), (x, y0)
    ]
    codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4,
             MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    ax.add_patch(PathPatch(MplPath(verts, codes), transform=ax.transAxes,
                           fill=False, lw=lw, color=color, clip_on=False))


def arrow(ax, xy0, xy1, color="black", lw=1.3, ms=14):
    ax.add_patch(FancyArrowPatch(xy0, xy1, transform=ax.transAxes,
                                 arrowstyle="->", mutation_scale=ms,
                                 lw=lw, color=color, clip_on=False))


def line(ax, x0, x1, y, **kw):
    ax.plot([x0, x1], [y, y], transform=ax.transAxes, clip_on=False, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/figures.yml")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
    fc = cfg["colors"]
    f1 = cfg["fig1"]
    W, H = f1["width_in"], f1["height_in"]
    fs_title = f1.get("title_font_size", 17)
    fs_time = f1.get("time_font_size", 15)
    fs_section = f1.get("section_font_size", 13)
    fs_body = f1.get("body_font_size", 10.5)
    fs_bottom = f1.get("bottom_heading_font_size", 10.4)

    fig, ax = plt.subplots(figsize=(W, H))
    ax.set_axis_off()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # Top estimand and timeline.
    ax.text(.50, .975, r"$S_{i,t_0}\ \longrightarrow\ A_{i,t_1}\ \longrightarrow\ Y_{i,t_2}$",
            ha="center", va="top", fontsize=fs_title, family="serif")
    ax.text(.50, .925, r"$t_0<t_1<t_2$", ha="center", va="top", fontsize=fs_time, family="serif")
    line(ax, .21, .94, .845, color="black", lw=1.2)
    arrow(ax, (.93,.845), (.955,.845), lw=1.2)
    xs=[.275,.55,.805]
    for x,t in zip(xs,[r"$t_0$",r"$t_1$",r"$t_2$"]):
        ax.scatter([x],[.845], s=38, color="black", transform=ax.transAxes, zorder=3)
        ax.plot([x,x],[.795,.845], ls=(0,(5,4)), lw=1, color="black", transform=ax.transAxes)
        ax.text(x, .862, t, ha="center", va="bottom", fontsize=13, family="serif")

    # Top left firm information.
    ax.text(.025,.775,r"Firm information $Z_i$",fontsize=13.5,fontweight="bold",family="serif")
    left_items=["Firm characteristics",r"Operating conditions at $t_0$","Survey structure"]
    for j,txt in enumerate(left_items):
        y=.735-j*.052
        ax.text(.025,y,u"•",fontsize=14,va="center",family="serif")
        ax.text(.04,y,txt,fontsize=11.5,va="center",family="serif")
    brace(ax,.171,.653,.748,"right")
    arrow(ax,(.184,.700),(.222,.700),lw=1.2)

    # Shock block.
    ax.text(.243,.748,r"Initial sales shock  $S_i$",fontsize=fs_section,color=fc["shock"],fontweight="bold",family="serif")
    ax.text(.243,.703,"Sales decline  (continuous intensity)",fontsize=fs_body,family="serif")
    ax.plot([.235,.405],[.655,.655],color=fc["shock"],lw=1.4,transform=ax.transAxes)
    arrow(ax,(.393,.655),(.414,.655),color=fc["shock"],lw=1.4,ms=12)
    vals=[0,25,50,75,80]; poss=np.linspace(.235,.405,5)
    for x,v in zip(poss,vals):
        ax.plot([x,x],[.646,.664],color=fc["shock"],lw=1,transform=ax.transAxes)
        ax.text(x,.632,f"{v}%",ha="center",va="top",fontsize=9.5,family="serif")
    arrow(ax,(.435,.675),(.482,.675),lw=1.2)

    # Treatment block.
    ax.text(.505,.748,r"Digital-channel expansion  $A_i$",fontsize=fs_section,color=fc["treatment"],fontweight="bold",family="serif")
    brace(ax,.505,.605,.698,"left")
    ax.text(.535,.687,r"$A_i=1$  (expanded)",fontsize=11.5,color=fc["treatment"],family="serif")
    ax.text(.535,.625,r"$A_i=0$  (not expanded)",fontsize=11.5,color=fc["treatment"],family="serif")

    # Outcomes.
    ax.text(.755,.748,r"Recovery outcomes  $Y_i$",fontsize=fs_section,color=fc["treatment"],fontweight="bold",family="serif")
    brace(ax,.755,.555,.695,"right")
    out=["Sales recovery","Continued operation","Employment retention","Input recovery"]
    for j,txt in enumerate(out):
        y=.685-j*.04
        ax.text(.765,y,u"•",fontsize=13,va="center",family="serif")
        ax.text(.78,y,txt,fontsize=10.7,va="center",family="serif")

    ax.text(.50,.535,
            "Include only firms with no prior reported digital-channel expansion; retain the earliest valid three-wave window for each firm",
            ha="center",fontsize=10.3,family="serif")
    line(ax,.015,.985,.512,color="black",lw=1.0,ls=(0,(5,3)))

    # Bottom: panel observations.
    ax.text(.03,.45,"Panel observations",fontsize=11.2,fontweight="bold",family="serif")
    ax.add_patch(FancyBboxPatch((.025,.15),.12,.26,boxstyle="round,pad=0.012",
                                transform=ax.transAxes,fill=False,lw=1.0))
    ax.text(.085,.38,r"$(S_i,Z_i,A_i,Y_i),$",ha="center",fontsize=10.3,family="serif")
    ax.text(.085,.35,r"$i=1,\ldots,N$",ha="center",fontsize=10.3,family="serif")
    x0=.055; y0=.17; cw=.026; ch=.04
    for r in range(3):
        for c in range(3):
            ax.add_patch(Rectangle((x0+c*cw,y0+r*ch),cw,ch,transform=ax.transAxes,fill=False,lw=.55))
    ax.text(x0+1.5*cw,y0+3*ch+.025,r"$t_0\quad t_1\quad t_2$",ha="center",fontsize=8.5,family="serif")
    ax.text(.03,.275,r"$i=1$",fontsize=8.5,family="serif")
    ax.text(.03,.235,r"$i=2$",fontsize=8.5,family="serif")
    ax.text(.035,.195,r"$\vdots$",fontsize=9,family="serif")
    ax.text(.03,.16,r"$i=N$",fontsize=8.5,family="serif")
    arrow(ax,(.15,.285),(.18,.285),lw=1.1)

    # Cross-fitting matrix.
    ax.text(.195,.45,r"Cross-fitting  $(K=5)$",fontsize=fs_bottom,fontweight="bold",family="serif")
    mx0=.215; my0=.245; cs=.02
    for r in range(5):
        for c in range(5):
            rect=Rectangle((mx0+c*cs,my0+(4-r)*cs),cs,cs,transform=ax.transAxes,
                           facecolor="white",edgecolor="black",lw=.55)
            ax.add_patch(rect)
            if r==c:
                rect.set_facecolor("#eef4ff")
                for k in np.linspace(0,cs*1.6,5):
                    ax.plot([mx0+c*cs+k-cs*.4,mx0+c*cs+k],
                            [my0+(4-r)*cs,my0+(4-r)*cs+cs],
                            transform=ax.transAxes,color=fc["treatment"],lw=.45,clip_on=True)
    ax.text(.205,.412,"Fold",fontsize=9.5,family="serif")
    for c in range(5): ax.text(mx0+(c+.5)*cs,.412,str(c+1),ha="center",fontsize=8,family="serif")
    for r in range(5): ax.text(mx0-.012,my0+(4-r+.35)*cs,str(r+1),ha="right",fontsize=8,family="serif")
    ax.text(.255,.185,"Out-of-fold prediction",ha="center",fontsize=9.8,family="serif")
    arrow(ax,(.32,.285),(.35,.285),lw=1.1)

    # DR pseudo outcome mini scatter.
    ax.text(.405,.45,"Doubly robust\northogonal signal $\\psi_i$",ha="center",fontsize=10.8,fontweight="bold",family="serif")
    rng=np.random.default_rng(2)
    sx=np.linspace(.36,.455,44); sy=.285+rng.normal(0,.025,len(sx))
    ax.scatter(sx,sy,s=8,color=fc["treatment"],transform=ax.transAxes)
    line(ax,.35,.47,.285,color="black",lw=.8,ls=(0,(4,3)))
    arrow(ax,(.35,.18),(.35,.40),lw=1)
    arrow(ax,(.35,.18),(.47,.18),lw=1)
    ax.text(.425,.145,r"$S_i$",fontsize=10,family="serif")
    ax.text(.345,.404,r"$\psi_i$",fontsize=10,family="serif")
    arrow(ax,(.475,.285),(.505,.285),lw=1.1)

    # Conditional g(s,z) mini-curve.
    ax.text(.57,.45,r"Conditional effect  $g(s,z)$",ha="center",fontsize=10.8,fontweight="bold",family="serif")
    xx=np.linspace(.51,.655,100); yy=.27+.025*np.sin(np.linspace(0,3*np.pi,100))+.07*(xx-.51)/.145
    band=.03
    ax.fill_between(xx,yy-band,yy+band,color="#cfe1f7",alpha=.8,transform=ax.transAxes)
    ax.plot(xx,yy,color=fc["treatment"],lw=1.5,transform=ax.transAxes)
    for _ in range(45):
        x=rng.uniform(.51,.655); y=np.interp(x,xx,yy)+rng.normal(0,.025)
        ax.scatter([x],[y],s=6,color="#6ca8ed",alpha=.8,transform=ax.transAxes)
    line(ax,.51,.655,.285,color="black",lw=.7,ls=(0,(4,3)))
    arrow(ax,(.51,.18),(.51,.40),lw=1); arrow(ax,(.51,.18),(.66,.18),lw=1)
    ax.text(.57,.145,r"$s$  (initial sales shock intensity)",ha="center",fontsize=8.5,family="serif")
    ax.text(.505,.405,r"$g(s,z)$",fontsize=9.5,family="serif")
    arrow(ax,(.665,.285),(.69,.285),lw=1.1)

    # Standardization reference distribution.
    ax.text(.74,.45,"Standardized to the\ncommon-support\nreference distribution",ha="center",fontsize=9.7,family="serif")
    zx=rng.normal(.71,.012,35); zy=rng.uniform(.21,.36,35)
    ax.scatter(zx,zy,s=7,color="#777777",transform=ax.transAxes)
    ax.text(.735,.25,r"$\int$",fontsize=35,family="serif")
    ax.text(.75,.22,r"$dF_{\mathcal{C}}(z)$",fontsize=10,family="serif")
    ax.text(.705,.16,r"$Z$",fontsize=10,family="serif")
    arrow(ax,(.78,.285),(.81,.285),lw=1.1)

    # Final tau(s) curve.
    ax.text(.89,.45,r"Recovery effect  $\tau(s)$",ha="center",fontsize=10.8,fontweight="bold",family="serif")
    xx=np.linspace(.835,.965,100); yy=.275+.018*np.sin(np.linspace(0,3*np.pi,100))+.035*(xx-.835)/.13
    ax.fill_between(xx,yy-.03,yy+.03,color="#dceff1",alpha=.9,transform=ax.transAxes)
    ax.plot(xx,yy,color=fc["teal"],lw=1.6,transform=ax.transAxes)
    line(ax,.83,.968,.285,color="black",lw=.7,ls=(0,(4,3)))
    arrow(ax,(.83,.18),(.83,.40),lw=1); arrow(ax,(.83,.18),(.97,.18),lw=1)
    ax.text(.90,.145,r"$s$  (initial sales shock intensity)",ha="center",fontsize=8.5,family="serif")
    ax.text(.835,.405,r"$\tau(s)$",fontsize=9.5,family="serif")

    ax.text(.50,.075,r"$\tau(s)=\int g(s,z)\,dF_{\mathcal{C}}(z)$",
            ha="center",fontsize=19,color=fc["treatment"],family="serif")

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight", pad_inches=.08)
    plt.close(fig)

if __name__ == "__main__":
    main()
