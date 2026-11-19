-- Gives every numbered rule of a specification an anchor of its own, and turns
-- every mention of a rule, in either specification or rationale, into a link
-- to it. Rules starting with S are the data file format's; rules starting
-- with A or T are the analysis format's. Links to files that are not pages of
-- this site go to the repository.
local REPO = "https://github.com/Katsevich-Lab/sceptre-hub/blob/main/"
local FOLDER = { S = "perturb-seq-data-format", A = "sceptre-analysis-format", T = "sceptre-analysis-format" }
local FOLDER_SET = { ["perturb-seq-data-format"] = true, ["sceptre-analysis-format"] = true }
local PAGES = { ["SPEC.md"] = true, ["RATIONALE.md"] = true, ["overview.qmd"] = true }

local input = quarto.doc.input_file:gsub("\\", "/")
local here = input:match("([^/]+)/[^/]+$")          -- the folder of the page being rendered
local is_spec = input:match("SPEC%.md$") ~= nil
local is_overview = input:match("overview%.qmd$") ~= nil

local function rule_url(id)
  local folder = FOLDER[id:sub(1, 1)]
  local page = (folder == here) and "SPEC.html" or ((FOLDER_SET[here] and "../" or "") .. folder .. "/SPEC.html")
  if is_spec and folder == here then page = "" end
  return page .. "#" .. id
end

-- text such as "(S8.3)," becomes "(", a link, ")," ; a group such as "T5" is left alone
local function link_rules(inline)
  local before, id, after = inline.text:match("^(.-)([AST]%d+%.%d+)(.*)$")
  if not id or before:match("[%w_/]$") or after:match("^[%w_]") then return nil end
  return { pandoc.Str(before), pandoc.Link(id, rule_url(id), "", { class = "rule-ref" }), pandoc.Str(after) }
end

-- the first cell of a rule's row: the anchor the links point to
local function anchor_rows(table)
  if not is_spec then return nil end
  for _, body in ipairs(table.bodies) do
    for _, row in ipairs(body.body) do
      local cell = row.cells[1]
      local id = cell and pandoc.utils.stringify(cell.contents)
      if id and id:match("^[AST]%d+%.%d+$") then
        cell.contents = { pandoc.Plain({ pandoc.Link(id, "#" .. id, "", { id = id, class = "rule-id" }) }) }
      end
    end
  end
  return table
end

-- a relative link to something this site does not render goes to the repository
local function fix_link(link)
  -- a full address of one of this site's pages, as the overviews write them
  local folder, page, rest = link.target:match("^https://github%.com/Katsevich%-Lab/sceptre%-hub/blob/main/([%w%-]+)/(%u+)%.md(.*)$")
  if not folder then
    folder, rest = link.target:match("^https://raw%.githack%.com/Katsevich%-Lab/sceptre%-hub/main/([%w%-]+)/overview%.html(.*)$")
    page = folder and "overview"
  end
  if folder and FOLDER_SET[folder] then
    link.target = ((folder == here) and "" or ("../" .. folder .. "/")) .. page .. ".html" .. rest
    return link
  end
  local target = link.target
  if target:match("^%a+:") or target:match("^#") then return nil end
  local path = target:gsub("#.*$", "")
  if PAGES[path:match("[^/]*$")] then return nil end
  local parts = {}
  for part in (here .. "/" .. path):gmatch("[^/]+") do
    if part == ".." then table.remove(parts) elseif part ~= "." then parts[#parts + 1] = part end
  end
  link.target = REPO:gsub("/blob/", path:match("/$") and "/tree/" or "/blob/") .. table.concat(parts, "/")
  return link
end

-- the browser tab and the search results show the page's own first heading
local function title_from_heading(doc)
  if doc.meta.title then return nil end
  for _, block in ipairs(doc.blocks) do
    if block.t == "Header" and block.level == 1 then
      doc.meta.pagetitle = pandoc.utils.stringify(block.content)
      return doc
    end
  end
end

if not FOLDER_SET[here] then return { { Str = link_rules } } end   -- the home page
return { { Pandoc = title_from_heading }, { Table = anchor_rows }, { Link = fix_link }, { Str = link_rules } }
