def preamble():
    preamble_text = r"""\usepackage[utf8]{inputenc}
\usepackage{natbib}
\usepackage{graphicx}
\usepackage[margin=1.0in]{geometry}
\usepackage{enumerate}
\usepackage[english]{babel}
\usepackage [autostyle, english = american]{csquotes}
\usepackage{amsmath}
\usepackage{siunitx}
\usepackage{setspace}
\usepackage{empheq}
\usepackage{textcomp}
\usepackage[labelfont=bf]{caption}
\usepackage{tikz}
\usepackage[superscript,biblabel]{cite}
\usepackage{hyperref}
\usepackage{achemso}
\usepackage{gensymb}
\usepackage{upgreek}
\urlstyle{same}
\graphicspath{{Images/}}
\DeclareSIUnit{\calorie}{cal}
\DeclareSIUnit{\kcal}{\kilo\calorie}
\sisetup{per-mode=symbol,retain-explicit-plus}
\newcommand{\beginsupplement}{%
        \setcounter{table}{0}
        \renewcommand{\thetable}{S\arabic{table}}%
        \setcounter{figure}{0}
        \renewcommand{\thefigure}{S\arabic{figure}}%
     }
"""
    return preamble_text