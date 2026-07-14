#pragma once
#include "matplot/matplot.h"

static void makePlot(int range)
{
    using namespace matplot;
    plot(   );
    xlim(0,range);
    ylim(0,range);
    xlabel("Distance from 3' end of  sequence read");
    ylabel("Mismatch frequency");
}
