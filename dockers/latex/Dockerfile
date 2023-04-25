FROM "luck/latex-base:latest"
#
# I know "scheme-full" is tempting, but three hours later you may have wished you picked out the packages...
#
# 1. Try using
#
#  tlmgr search --global --file <missing thing>
#
# to figure out what you need to install.
#
# 2. append extra RUN lines so docker can cache your build steps instead of re-running all of them.
#
# 3. Once you have the packages, combine them into one tlmgr install in sorted order to help docker
#    cache them once and make it easier on your collaborators...
#
# 4. You can nullify multiline RUN lines by replacing inserting an 'echo skip' at the begining:
#
#    # enabled
#    RUN tlmgr install \
#        <package1> \
#        <package2>
#
#    # disabled
#    RUN echo skip tlmgr install \
#        <package1> \
#        <package2>
#
#    This format makes it easier to manage a package list in sorted order, which is the
#    recommendation after you figure out all the packages to install through multiple
#    RUN commands.
#

#
# scheme-full (5+GB, 2+ hours) - comment out to disable
#
RUN echo skip tlmgr install scheme-full

#
# min (remove echo skip below and comment out scheme-full above to use)
#
RUN echo skip tlmgr install \
    collection-fontsrecommended

#
# tufte (remove echo skip below and comment out scheme-full above to use)
#
RUN tlmgr install \
    asymptote \
    bera \
    booktabs \
    catchfile \
    ccicons \
    collection-fontsrecommended \
    fancyvrb \
    fpl \
    hardwrap \
    l3kernel \
    l3packages \
    lipsum \
    listings \
    mathpazo \
    microtype \
    ms \
    nth \
    paralist \
    palatino \
    placeins \
    ragged2e \
    relenc \
    setspace \
    siunitx \
    soul \
    textcase \
    titlesec \
    tufte-latex \
    units \
    xcolor \
    xifthen \
    xindy \
    xkeyval

RUN tlmgr install \
    etoolbox

# Install sympy

RUN apt update
RUN apt install -y python3-pip
RUN pip3 install mpmath sympy

RUN rm -rf /usr/local/texlive/2023/texmf-dist/tex/latex/units
