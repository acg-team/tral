:orphan:

Add a class inheritance diagram within .rst files:

::
.. inheritance-diagram:: tral.repeat.repeat tral.sequence.sequence tral.repeat_list.repeat_list tral.hmm.hmm
   :parts: 2


Create a UML-like diagram from the codebase with pylint's pyreverse:
::

    pyreverse -o svg -p Tral tral/


Prettify the produces graphs following http://www.hokstad.com/making-graphviz-output-pretty-with-xsl
::

    xsltproc notugly.xsl _static/packages_Tral.svg > _static/packages_Tral_notugly.svg
    xsltproc notugly.xsl _static/classes_Tral.svg > _static/classes_Tral_notugly.svg