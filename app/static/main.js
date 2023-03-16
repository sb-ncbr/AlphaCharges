"use strict";

let molstar;

function init_results(structure_url, id) {
    (async () => {
        molstar = await MolstarPartialCharges.create("root");
        await load(structure_url, id);
    })().then(
        () => {},
        (error) => {
            console.error("Mol* initialization ❌", error);
        }
    );
}

function init_wrong_structure(structure_url, problematicAtoms) {
    const parsedAtoms = JSON.parse(problematicAtoms);
    addProblematicAtoms(parsedAtoms);

    (async () => {
        molstar = await MolstarPartialCharges.create("root");
        await molstar.load(structure_url, "pdb");
        await molstar.type.ballAndStick();
        await molstar.color.default();
    })().then(
        () => {},
        (error) => {
            console.error("Mol* initialization ❌", error);
        }
    );
}

async function load(structure_url, id) {
    const first_example = "P34712_7.2_4";

    await molstar.load(structure_url);

    if (id === first_example) {
        document.getElementById("view_surface").setAttribute("checked", "true");
        await molstar.type.surface();
    } else {
        document
            .getElementById("colors_relative")
            .setAttribute("checked", "true");
        await molstar.type.default();
    }

    updateRelativeColor();
    mountTypeControls();
    mountColorControls();
}

function mountTypeControls() {
    const cartoon = document.getElementById("view_cartoon");
    const surface = document.getElementById("view_surface");
    const bas = document.getElementById("view_bas");
    if (!cartoon || !surface || !bas) return;
    cartoon.onclick = async () => await molstar.type.default();
    surface.onclick = async () => await molstar.type.surface();
    bas.onclick = async () => await molstar.type.ballAndStick();
}

function mountColorControls() {
    const structure = document.getElementById("colors_structure");
    const alphafold = document.getElementById("colors_alphafold");
    const relative = document.getElementById("colors_relative");
    const absolute = document.getElementById("colors_absolute");
    const range = document.getElementById("max_value");
    if (!structure || !relative || !absolute || !alphafold) return;
    structure.onclick = async () => await updateDefaultColor();
    alphafold.onclick = async () => await updateAlphaFoldColor();
    relative.onclick = async () => await updateRelativeColor();
    absolute.onclick = async () => await updateAbsoluteColor();
    range.oninput = async () => await updateRange();
}

async function updateDefaultColor() {
    const input = document.getElementById("max_value");
    if (!input) return;
    input.setAttribute("disabled", "true");
    await molstar.color.default();
}

async function updateAlphaFoldColor() {
    const input = document.getElementById("max_value");
    if (!input) return;
    input.setAttribute("disabled", "true");
    await molstar.color.alphaFold();
}

async function updateRelativeColor() {
    const input = document.getElementById("max_value");
    if (!input) return;
    input.setAttribute("disabled", "true");
    const charge = await molstar.charges.getRelativeCharge();
    input.value = charge.toFixed(3);
    // ? updates viewer each time the radio button is clicked ?
    await molstar.color.relative();
}

async function updateAbsoluteColor() {
    const input = document.getElementById("max_value");
    if (!input) return;
    input.removeAttribute("disabled");
    await molstar.color.relative();
}

async function updateRange() {
    const input = document.getElementById("max_value");
    if (!input) return;
    const value = Number(input.value);
    if (isNaN(value)) return;
    await molstar.color.absolute(value);
}

function addProblematicAtoms(problematicAtoms) {
    const span = document.getElementById("problematic_atoms");
    if (!span) return;
    Object.keys(problematicAtoms).forEach((id, i) => {
        const button = createProblematicAtomButton(
            id,
            problematicAtoms[id].key
        );
        const tooltip = createProblematicAtomTooltip(
            problematicAtoms[id].message
        );
        span.appendChild(button);
        span.appendChild(tooltip);
        if (i < Object.keys(problematicAtoms).length - 1) {
            const text = document.createTextNode(", ");
            span.appendChild(text);
        }
    });
}

function createProblematicAtomButton(id, key) {
    const button = document.createElement("button");
    button.id = id;
    button.className = "btn btn-link p-0 font-weight-bold";
    button.onclick = () => molstar.visual.focus(key);
    button.textContent = id;
    return button;
}

function createProblematicAtomTooltip(message) {
    const tooltip = document.createElement("i");
    tooltip.className = "bi bi-question";
    tooltip.setAttribute("data-toggle", "tooltip");
    tooltip.setAttribute("data-placement", "top");
    tooltip.setAttribute("title", message);
    return tooltip;
}
