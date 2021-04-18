import React from "react";

const styles = {
  borderBottom: "2px solid #eee",
  background: "#fafafa",
  margin: ".75rem auto",
  padding: ".6rem 1rem",
  maxWidth: "500px",
  borderRadius: "7px"
};

export default ({
  data: {
    recipeName: name,
    recipeIngredients: ingredients,
    recipeMethods: methods
  },
  handleClick,
  handleDelete
}) => {
  return (
    <div style={styles}>
      <h2>{name}</h2>
      <p>{ingredients}</p>
      <p>{methods}</p>

      <div className="box-footer">
        <button
          style={{ marginRight: 10 }}
          className="btn btn-primary"
          type="button"
          onClick={handleClick}
        >
          update
        </button>
        <button className="btn btn-danger" type="button" onClick={handleDelete}>
          Delete
        </button>
      </div>
    </div>
  );
};
