import React, { Component } from "react";
import { connect } from "react-redux";
import { updateRecipe, addRecipe } from "../../../../store/actions/recipes";
import { getPostById } from "../../../../services/recipes";

class NewRecipe extends Component {
  constructor() {
    super();
    this.state = {
      recipeName: "",
      recipeIngredients: "",
      recipeMethods: "",
      error: false
    };
  }

  async componentDidMount() {
    const { id } = this.props.match.params;
    if (id === "new-recipe") return;

    //fetch data of particular id
    try {
      const data = await getPostById(id);
      this.mapToState(data);
    } catch (error) {
      this.setState({ error });
    }
  }

  mapToState = ({ recipeIngredients, recipeMethods, recipeName }) => {
    this.setState({
      recipeIngredients,
      recipeMethods,
      recipeName
    });
  };

  handleChange = ({ currentTarget }) => {
    const { value, id } = currentTarget;
    this.setState({ [id]: value });
  };

  handleSubmit = async e => {
    e.preventDefault();
    const { recipeName, recipeIngredients, recipeMethods } = this.state;

    if (this.props.match.params.id === "new-recipe") {
      try {
        await this.props.addRecipe({
          recipeName,
          recipeIngredients,
          recipeMethods
        });
        this.props.history.push("/mainpage/admin/recipe");
      } catch (error) {
        this.setState({ error });
      }
    } else {
      try {
        this.props.updateRecipe(
          {
            recipeName,
            recipeIngredients,
            recipeMethods
          },
          this.props.match.params.id
        );

        this.props.history.push("/mainpage/admin/recipe");
      } catch (error) {
        this.setState({ error });
      }
    }
  };

  render() {
    const { recipeName, recipeIngredients, recipeMethods, error } = this.state;
    const { id } = this.props.match.params;
    // const { error } = this.props;

    if (error) return <h1>{error.message}</h1>;
    return (
      <div className="col-sm-6">
        <div className="box box-info">
          <div className="box-header with-border">
            <h3 className="box-title">Recipe Form</h3>
          </div>
          <form className="form-horizontal" onSubmit={this.handleSubmit}>
            <div className="box-body">
              <div className="form-group">
                <label htmlFor="recipeName" className="col-sm-2 control-label">
                  Name
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter recipe name"
                    type="text"
                    className="form-control"
                    id="recipeName"
                    defaultValue={recipeName}
                    onChange={this.handleChange}
                  />
                </div>
              </div>

              <div className="form-group">
                <label
                  htmlFor="recipeIngredients"
                  className="col-sm-2 control-label"
                >
                  Ingredients
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter recipe ingredients"
                    type="text"
                    className="form-control"
                    id="recipeIngredients"
                    defaultValue={recipeIngredients}
                    onChange={this.handleChange}
                  />
                </div>
              </div>

              <div className="form-group">
                <label
                  htmlFor="recipeMethods"
                  className="col-sm-2 control-label"
                >
                  Method
                </label>
                <div className="col-sm-10">
                  <input
                    placeholder="enter recipe method"
                    rows="5"
                    type="text"
                    className="form-control"
                    id="recipeMethods"
                    defaultValue={recipeMethods}
                    onChange={this.handleChange}
                  />
                </div>
              </div>
            </div>
            <div className="box-footer">
              <button className="btn btn-primary pull-right">
                {id === "new-recipe" ? "Submit" : "Update"}
              </button>
            </div>
          </form>
          {error && <div>{error}</div>}
        </div>
      </div>
    );
  }
}

const mapStateToProps = ({ recipeManagement }) => ({
  error: recipeManagement.recipeList.error
});

const mapDsipatchToProps = dispatch => ({
  addRecipe: data => dispatch(addRecipe(data)),
  updateRecipe: (data, id) => dispatch(updateRecipe(data, id))
});

export default connect(
  mapStateToProps,
  mapDsipatchToProps
)(NewRecipe);
