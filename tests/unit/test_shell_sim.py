"""The SampleEditor that backs both ``.lcm`` parsing and interactive SIM.

The behaviour under test is stated in sim.htm: ``LAYER n`` *navigates* rather
than appending, there is always one blank layer at the bottom, and a layer left
at zero thickness "disappears automatically if the active layer changes".
"""

from __future__ import annotations

import pytest

from pyrump.script.lcm import SampleEditor, parse_lcm, write_lcm

ITO = """Sim Reset
Layer 1
 Thick 151 ITO
 Composition In 2 O 3 Sn 0.1 /
Next
 Thick 10 um
 Composition O 4 C 14 H 10 /
Maxpth 200
"""


def run(*lines: str) -> SampleEditor:
    editor = SampleEditor()
    for line in lines:
        editor.execute(line)
    return editor


def test_parse_lcm_is_a_loop_over_the_editor():
    script = parse_lcm(ITO)
    assert len(script.layers) == 2
    assert script.layers[0].composition == {"In": 2.0, "O": 3.0, "Sn": 0.1}
    assert script.layers[0].unit == "ITO"
    assert script.layers[1].unit == "um"
    assert script.maxpth == 200.0


def test_layer_navigates_rather_than_appending():
    editor = run(
        "Layer 1", "Thick 100 A", "Composition Si 1 /",
        "Next", "Thick 200 A", "Composition Au 1 /",
        # Going back to layer 1 must not create a third layer.
        "Layer 1", "Thick 150 A",
    )
    script = editor.finish()
    assert len(script.layers) == 2
    assert script.layers[0].thickness == 150.0
    assert script.layers[1].thickness == 200.0


def test_next_at_the_bottom_lands_on_the_blank_layer():
    editor = run("Layer 1", "Thick 100 A", "Composition Si 1 /", "Next")
    # The pointer is on the blank layer; nothing is committed until it is given
    # substance, so finishing here yields one layer.
    assert len(editor.finish().layers) == 1


def test_a_zero_thickness_layer_disappears_when_the_pointer_moves():
    editor = run(
        "Layer 1", "Thick 100 A", "Composition Si 1 /",
        "Next",           # blank layer
        "Layer 1",        # moving away prunes it
    )
    assert len(editor.finish().layers) == 1


def test_open_inserts_above_the_current_layer():
    editor = run(
        "Layer 1", "Thick 100 A", "Composition Si 1 /",
        "Next", "Thick 200 A", "Composition Au 1 /",
        "Layer 2", "Open", "Thick 50 A", "Composition C 1 /",
    )
    script = editor.finish()
    assert [layer.thickness for layer in script.layers] == [100.0, 50.0, 200.0]
    assert list(script.layers[1].composition) == ["C"]


def test_a_layer_keeps_its_composition_when_only_thickness_is_zero():
    # Composition alone is enough to keep a layer, matching parse_lcm's original
    # "thickness > 0 or composition" rule.
    editor = run("Layer 1", "Composition Si 1 /", "Next", "Thick 10 A")
    assert len(editor.finish().layers) == 2


def test_reset_clears_the_sample():
    editor = run("Layer 1", "Thick 100 A", "Composition Si 1 /", "Sim Reset")
    assert editor.finish().layers == []


def test_unknown_commands_are_collected_not_raised():
    editor = run("Layer 1", "Thick 100 A", "Composition Si 1 /", "Plot 1", "Wiggle 3")
    assert editor.finish().ignored == ["Plot 1", "Wiggle 3"]


def test_an_unknown_equation_raises():
    with pytest.raises(ValueError, match="unknown equation"):
        run("Layer 1", "Thick 100 A", "Equation nonsense 1 2 3")


def test_editing_round_trips_through_write_lcm():
    script = run(
        "Sim Reset",
        "Layer 1", "Thick 151 ITO", "Composition In 2 O 3 Sn 0.1 /",
        "Next", "Thick 10 um", "Composition O 4 C 14 H 10 /",
        "Maxpth 200",
    ).finish()
    assert parse_lcm(write_lcm(script)).layers[0].composition == (
        script.layers[0].composition
    )
