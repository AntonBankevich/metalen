Error = ""
Ready = False
try:
    import matplotlib
    if matplotlib.__version__.startswith('0') or matplotlib.__version__.startswith('1.0'):
        Error = "Matplotlib version is old. Please use version 1.1 or higher."
        Ready = False
    else:
        import matplotlib.pyplot as plt
        import matplotlib.backends.backend_pdf
        Ready = True
except Exception:
    Error = "Matplotlib version is not istalled or corrupted."
    Ready = False

import math

import sys


def DrawHeightHistogram(ax, heights, color = "black"):
    curvex, curvey = HeightsToCurve(heights)
    #plt.plot(x, y, cols[0])
    ax.plot(curvex, curvey, color)
    # map(lambda (d, c): ax.plot(d['X'], d['Y'], c, markersize=4, markeredgewidth=0.0), zip(df, cols)[0:])
    ax.set_yscale("log")
    scale = ""
    # if df[0]["m"] == 1000000:
    #     scale = "(Mb)"
    # elif df[0]["m"] == 1000000000:
    #     scale = "(Gb)"
    ax.set_xlabel("Cumulative length" + scale)
    ax.set_ylabel("Height", color = color)
    ax.set_xlim(min(curvex) - 1, max(curvex)*1.1)
    ax.tick_params(axis='y', colors=color)
    # plt.xticks(list(plt.xticks()[0]) + [curvex[-1]])
    mi = min(curvey)
    mi1 = math.pow(10, int(math.log10(mi)))
    if mi1 / mi > 1.5:
        mi = mi1 / 10
    ma = max(curvey)
    ma = math.pow(10, int(math.log10(ma)))
    ax.set_ylim(mi, ma)
    # ax.set_title(name)

def DrawFigure(fname, heights):
    ax = plt.subplot()
    DrawHeightHistogram(ax, heights)
    plt.tight_layout()
    pp = matplotlib.backends.backend_pdf.PdfPages(fname)
    pp.savefig()
    pp.close()
    plt.clf()



def HeightsToCurve(heights):
    heights = sorted(heights)[::-1]
    cur = 0
    all_points = []
    for h in heights:
        cur += 1 / h
        all_points.append((cur / len(heights), h))
    sum = all_points[-1][0]
    step = max(1, len(all_points) / 2000)
    cnt = 0
    prev = 0
    resx = []
    resy = []
    for x, y in all_points:
        if cnt % step == 0 or x - prev > sum / 2000:
            resx.append(x)
            resy.append(y)
            prev = x
        cnt += 1
    return resx, resy

if __name__ == "__main__":
    heights = map(lambda x: float(x.split()[-1]), open(sys.argv[1], "r").readlines())
    DrawFigure("fig.pdf", heights)