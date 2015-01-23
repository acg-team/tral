.. _significance_test:

Perform statistical significance test of tandem repeats.
========================================================

Here you learn run statistical significance tests on your candidate set of tandem repeats.
TRAL currently provides a variety of model-based and none-model based tests. They are
detailed in a :ref:`NAR publication (2012) <publications>`.

[INCOMPLETE]

Read in tandem repeat annotations.
---------------------------------------

::

    import os, pickle
    from tandemrepeats import repeat_list

    path_to_tandemrepeats = "/my/path/to/TandemRepeats/tandemrepeats"
    fRepeat_Pickle = os.path.join(path_to_tandemrepeats,"test","HIV-1_388796.pickle")

    with open(fRepeat_Pickle, "rb") as fh:
        dTandem_repeats = pickle.load(fh)

    print(dTandem_repeats)



Perform a statistical significance test
---------------------------------------

Calculate the statistical score, and the corresponding p-Value with parameters defined in
the :ref:`configuration file <configure>`:
::

    for iRepeatList in dTandem_repeats.values():
        if iRepeatList:
            for iTandemRepeat in iRepeatList.repeats:
                iTandemRepeat.calculate_pValues()

