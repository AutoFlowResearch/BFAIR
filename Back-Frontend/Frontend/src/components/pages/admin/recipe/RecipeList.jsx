import React, { Component } from 'react';
import { connect } from 'react-redux';
import { getAllRecipes, deleteRecipe } from '../../../../store/actions/recipes';
import Card from '../../../common/Card';

class RecipeList extends Component {
  componentDidMount() {
    try {
      this.props.getAllRecipes();
    } catch (error) {}
  }

  handleNewRecipe = () => {
    this.props.history.push('/mainpage/admin/recipe/new-recipe');
  };

  handleRecipe = (id) => {
    this.props.history.push(`/mainpage/admin/recipe/${id}`);
  };

  handleDelete = (id) => {
    this.props.deleteRecipe(id);
  };

  render() {
    const { data, error, isFetching } = this.props;
    if (error) return <h3>{error.message}</h3>;
    if (isFetching) return <h3>Loading</h3>;

    return (
      <div>
        <button
          className="btn btn-primary"
          style={{ margin: 9 }}
          onClick={this.handleNewRecipe}
        >
          chandu Recipe
        </button>
        {data.map((recipe) => (
          <Card
            key={recipe._id}
            handleClick={() => this.handleRecipe(recipe._id)}
            handleDelete={() => this.handleDelete(recipe._id)}
            data={recipe}
          />
        ))}
      </div>
    );
  }
}

const mapStateToProps = (state) => ({
  data: state.recipeManagement.recipeList.data,
  error: state.recipeManagement.recipeList.error,
  isFetching: state.recipeManagement.recipeList.isFetching,
});

const mapDispatchToProps = (dispatch) => ({
  getAllRecipes: () => dispatch(getAllRecipes()),
  deleteRecipe: (id) => dispatch(deleteRecipe(id)),
});

export default connect(mapStateToProps, mapDispatchToProps)(RecipeList);
